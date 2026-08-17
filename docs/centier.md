# 着丝粒鉴定 | CentIER Centromere Identification

一句话理解：**在 T2T 基因组上自动圈出每条染色体的「着丝粒」位置**——综合 k-mer 频率、串联重复、LTR 转座子等多重特征预测着丝粒区域，可选 Hi-C 数据进一步验证边界。

## 功能概述 | Overview { #overview }

- 封装 CentIER v3.0.1，用于 T2T（端粒到端粒）基因组的着丝粒识别与注释
- 多特征融合：k-mer 频率分析 + 串联重复 + LTR 反转座子
- 可选 GFF 注释提升精度；可选 Hi-C 矩阵验证边界
- 支持 Hi-C FASTQ 自动模式：内置 HiC-Pro 流程直接产出矩阵，再接 CentIER
- 输出着丝粒坐标、序列、单体、LTR 位置统计和 SVG 可视化图

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools centier -i genome.fa -o output_dir/
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| 着丝粒 | 染色体中间「收腰」的关键区域，细胞分裂时负责把染色体拉向两极 |
| T2T 基因组 | 「端粒到端粒」无缺口的完整组装，着丝粒这类重复区也被完整拼出 |
| 单体（monomer） | 着丝粒里反复出现的重复单元，像同一句口号重复成千上万遍 |
| 串联重复 | 一段序列紧挨着自己重复多次，着丝粒区尤其富集 |
| LTR 转座子 | 会「自我复制搬家」的重复元件，着丝粒区常富集特定类型（如 Ty3_gypsy） |
| k-mer 频率 | 把序列切成 k 长度小片段后统计多样性；着丝粒区重复度高、k-mer 种类骤降 |
| Hi-C | 测「染色体在三维空间里谁和谁挨得近」的技术，着丝粒区相互作用信号会减弱 |

## 输入 | Input { #input }

- `-i` 输入：基因组 FASTA（必需）
- 可选注释：`--gff`（GFF/GTF）
- Hi-C 两种给法（二选一）：
  - 现成矩阵：`--matrix1` + `--matrix2` + `--bed1` + `--bed2`（四个必须一起给，100kb 和 20kb 分辨率各一对）
  - 原始 FASTQ 自动模式：`-1 fastq_r1 -2 fastq_r2`（提供即自动跑 HiC-Pro 产出矩阵）

染色体命名建议符合 `ChrN` 格式（如 `Chr1`、`So120_Chr1`）；`Chr_1` / `scaffold_1` 会让 CentIER 崩溃，默认只警告、加 `--strict-chrname` 则直接中止。

## 参数说明 | Parameters { #parameters }

### 软件与注释 | Software and annotation

**通俗理解|In plain words:** `--centier-path` 是 CentIER 安装目录（默认 `~/software/CentIER/CentIER-main`），按默认装就不用改；`--gff` 是可选的注释文件，能提高预测精度但非必需。**没有注释文件时 CentIER 也能正常跑。**

相关参数：`--centier-path`、`--gff`。

### Hi-C 数据 | Hi-C data

**通俗理解|In plain words:** 现成矩阵模式需要四个文件齐活（100kb 矩阵+bed、20kb 矩阵+bed）；FASTQ 自动模式（`-1/-2`）会自动跑 HiC-Pro 产出矩阵。`--hic-matrix-type` 选用 raw 还是 iced 归一化矩阵，`--force-hicpro` 强制重跑 HiC-Pro。**没有 Hi-C 数据也能跑，只是少了边界验证这一环。**

相关参数：`--matrix1`、`--matrix2`、`--bed1`、`--bed2`、`--fastq-r1`、`--fastq-r2`、`--genome-id`、`--restriction-enzyme`、`--bowtie2-idx`、`--bin-sizes`、`--max-memory`、`--force-hicpro`、`--hic-matrix-type`、`--strict-chrname`。

### 分析参数 | Analysis parameters

**通俗理解|In plain words:** `-k` 是 k-mer 长度；`-c` 中心容差决定「容许的中心位置偏差」；`--step-len` 是滑窗步长；`--signal-threshold` 和 `--mingap` 只在给 Hi-C 数据时用于信号判定。这些是 CentIER 原版参数，**默认值对大多数植物基因组够用，一般不用动**；染色体疑似多着丝粒（metapolycentric）时才加 `--mul-cents` 保留所有候选区。

相关参数：`--kmer-size`、`--center-tolerance`、`--step-len`、`--mul-cents`、`--mingap`、`--signal-threshold`。

### 运行控制 | Run control

**通俗理解|In plain words:** `--skip-dependency-check` 跳过开头的依赖检查（确认依赖齐全时用）；`--summary` 额外生成一份 JSON 摘要。`--step` 参数仅作兼容保留，实际仍会跑完整流程（见 FAQ）。

相关参数：`--skip-dependency-check`、`--summary`、`--step`。

## 分析流程 | Pipeline { #pipeline }

```text
基因组 FASTA
    │
    ├─ [Hi-C 自动模式] ChrN 命名预检 → HiC-Pro 产出矩阵 → 定位 100kb/20kb 矩阵
    │
    ▼
CentIER 鉴定（k-mer 频率 + 串联重复 + LTR 分类）
    │
    ├─ [有 Hi-C 数据] 用相互作用信号验证/精细化着丝粒边界
    │
    ▼
输出着丝粒坐标/序列/单体/LTR/SVG + 软件版本信息
```

## 输出 | Output { #output }

```text
输出目录/   （Hi-C 自动模式下为 02_centier/ 子目录）
├── {prefix}_centromere_range.txt         # 着丝粒区域坐标（染色体 起始 终止）
├── {prefix}_all_centromere_seq.txt       # 着丝粒区域序列
├── {prefix}_monomer_seq.txt              # 单体序列
├── {prefix}_monomer_in_centromere.txt    # 着丝粒内单体位置
├── {prefix}_ltr_position.txt             # LTR 位置
├── {prefix}_LTR_statistics.txt           # LTR 统计
├── {prefix}_draw_cen.svg                 # 着丝粒可视化图
├── centier_summary.json                  # 结果摘要（--summary 时）
├── 00_pipeline_info/software_versions.yml # 软件版本与参数
└── 99_logs/centier.log                   # 运行日志
```

Hi-C 自动模式额外产生 `01_hic_mapping/`（HiC-Pro 输出）。

## 结果解读 | Interpreting Results { #interpreting-results }

- `centromere_range.txt`：每行一个着丝粒（染色体 + 起止坐标），是核心结论，可直接用于后续分析或可视化
- `draw_cen.svg`：每条染色体上着丝粒区域的示意图，浏览器可打开；正常情况每条染色体应有一个清晰的主着丝粒区
- `LTR_statistics.txt`：各着丝粒区的 LTR 类型构成，植物着丝粒常富集 Ty3_gypsy 类
- 若有 Hi-C 数据：着丝粒区在矩阵上表现为信号减弱区，边界更可信
- 用 `--summary` 得到的 `centier_summary.json` 里 `centromere_count` 直接给出识别到的着丝粒数量

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 常规植物 T2T 基因组：`-i genome.fa -o out --summary`，其余默认
- 有现成 Hi-C 矩阵：补上 `--matrix1/2 --bed1/2` 增强边界可信度
- 只有 Hi-C FASTQ：`-1 r1.fq -2 r2.fq`（自动跑 HiC-Pro）
- 染色体名不是 ChrN 格式：先规范命名，或加 `--strict-chrname` 让程序显式中止而非静默出错
- 疑似多着丝粒物种：加 `--mul-cents`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 基因组FASTA文件｜Genome FASTA file path |
| `--output-dir, -o` | `./centier_output` |  | 输出目录｜Output directory |
| `--centier-path` | `~/software/CentIER/CentIER-main` |  | CentIER软件路径｜CentIER software path |
| `--gff` | — |  | GFF/GTF注释文件｜GFF/GTF annotation file (optional) |
| `--matrix1` | — |  | Hi-C矩阵文件(100000分辨率)｜Hi-C matrix file at 100000 resolution (optional) |
| `--matrix2` | — |  | Hi-C矩阵文件(200000分辨率)｜Hi-C matrix file at 200000 resolution (optional) |
| `--bed1` | — |  | Hi-C BED文件(对应matrix1)｜Hi-C BED file for matrix1 (optional) |
| `--bed2` | — |  | Hi-C BED文件(对应matrix2)｜Hi-C BED file for matrix2 (optional) |
| `--fastq-r1, -1` | — |  | Hi-C R1 FASTQ(启用自动模式)｜Hi-C R1 FASTQ (enables auto mode) |
| `--fastq-r2, -2` | — |  | Hi-C R2 FASTQ｜Hi-C R2 FASTQ |
| `--genome-id, -g` | — |  | 基因组ID｜Genome ID |
| `--restriction-enzyme` | `MboI` |  | 限制性内切酶｜Restriction enzyme |
| `--bowtie2-idx` | — |  | Bowtie2索引路径｜Bowtie2 index path |
| `--bin-sizes` | `100000 20000` |  | HiC-Pro bin大小｜HiC-Pro bin sizes |
| `--max-memory` | `200` | int | HiC-Pro最大内存GB｜HiC-Pro max memory GB |
| `--force-hicpro` | — |  | 强制重跑HiC-Pro｜Force rerun HiC-Pro |
| `--hic-matrix-type` | `raw` | raw/iced | Hi-C矩阵类型｜Hi-C matrix type |
| `--strict-chrname` | — |  | 染色体命名不符ChrN时中止｜Abort if chr naming not ChrN |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--kmer-size, -k` | `21` | int | K-mer大小｜K-mer size |
| `--center-tolerance, -c` | `15` | int | 中心容差｜Center tolerance |
| `--step-len` | `10000` | int | 步长｜Step length |
| `--mul-cents` | — |  | 保留所有潜在的着丝粒区域｜Retain all potential centromeric regions |
| `--mingap` | `2` | int | 最小Gap值｜Minimum gap value n*100000 |
| `--signal-threshold` | `0.7` | float | 信号阈值｜Signal threshold |
| `--step, -s` | — | IntRange | 运行指定步骤｜Run only specified step (1-6) |
| `--skip-dependency-check` | — |  | 跳过依赖检查｜Skip dependency check |
| `--summary` | — |  | 输出分析结果摘要｜Output analysis result summary |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file path |
| `-o, --output-dir` | `./centier_output` |  | 输出目录｜Output directory (default: ./centier_output) |
| `--centier-path` | `~/software/CentIER/CentIER-main` |  | CentIER软件路径｜CentIER software path (default: ~/software/CentIER/CentIER-main) |
| `--gff` | — |  | GFF/GTF注释文件｜GFF/GTF annotation file path (optional) |
| `--matrix1` | — |  | Hi-C矩阵文件(100000分辨率)｜Hi-C matrix file at 100000 resolution (optional) |
| `--matrix2` | — |  | Hi-C矩阵文件(200000分辨率)｜Hi-C matrix file at 200000 resolution (optional) |
| `--bed1` | — |  | Hi-C BED文件(对应matrix1)｜Hi-C BED file for matrix1 (optional) |
| `--bed2` | — |  | Hi-C BED文件(对应matrix2)｜Hi-C BED file for matrix2 (optional) |
| `-1, --fastq-r1` | — |  | Hi-C R1 FASTQ(提供即启用自动模式)｜Hi-C R1 FASTQ (enables auto mode) |
| `-2, --fastq-r2` | — |  | Hi-C R2 FASTQ｜Hi-C R2 FASTQ |
| `-g, --genome-id` | — |  | 基因组ID(bowtie2索引命名,默认从基因组文件名推导)｜Genome ID for bowtie2 index naming (default: derived from genome filename) |
| `--restriction-enzyme` | `MboI` |  | 限制性内切酶｜Restriction enzyme (default: MboI) |
| `--bowtie2-idx` | — |  | Bowtie2索引路径(默认自动建)｜Bowtie2 index path (auto-built if not given) |
| `--bin-sizes` | `100000 20000` |  | HiC-Pro bin大小(空格分隔)｜HiC-Pro bin sizes (default: 100000 20000) |
| `--max-memory` | `200` | int | HiC-Pro最大内存(GB)｜HiC-Pro max memory in GB (default: 200) |
| `--force-hicpro` | — | store_true | 强制重跑HiC-Pro｜Force rerun HiC-Pro |
| `--hic-matrix-type` | `raw` | raw/iced | Hi-C矩阵类型｜Hi-C matrix type (default: raw) |
| `--strict-chrname` | — | store_true | 染色体命名不符ChrN时中止｜Abort if chr naming not ChrN |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `-k, --kmer-size` | `21` | int | K-mer大小｜K-mer size (default: 21) |
| `-c, --center-tolerance` | `15` | int | 中心容差｜Center tolerance (default: 15) |
| `--step-len` | `10000` | int | 步长｜Step length (default: 10000) |
| `--mul-cents` | — | store_true | 保留所有潜在的着丝粒区域｜Retain all potential centromeric regions |
| `--mingap` | `2` | int | 最小Gap值｜Minimum gap value n*100000 (default: 2) |
| `--signal-threshold` | `0.7` | float | 信号阈值｜Signal threshold (default: 0.7) |
| `-s, --step` | — | 1/2/3/4/5/6 | 运行指定步骤｜Run only specified step (1-6) |
| `--skip-dependency-check` | — | store_true | 跳过依赖检查｜Skip dependency check |
| `--summary` | — | store_true | 输出分析结果摘要｜Output analysis result summary |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- CentIER v3.0.1（`centIER.py` + `bin/` 自带 hmmsearch、ltr_finder、trf409.linux64、REXdb.hmm、Ty3_gypsy.hmm）
- genometools（`gt`）与 LTR_retriever（须在 PATH 中，CentIER 必需）
- Python 包：pyfastx、numpy、pandas、scipy
- Hi-C 自动模式额外需要：HiC-Pro（默认 `~/software/HiC-Pro_v3.1.0`）、bowtie2-build

## 常见问题 | FAQ { #faq }

**Q1：报「未找到 centIER.py」或依赖工具缺失？**
确认 `--centier-path` 指向 CentIER 安装目录，且 `bin/` 下工具存在、`hmmsearch`/trf/ltr_finder 有可执行权限（`chmod +x`）；`gt` 和 LTR_retriever 需在 PATH（可 `conda install -n centier -c bioconda genometools`）。

**Q2：染色体命名报「not ChrN format」？**
把序列 ID 改成以 `ChrN` 结尾（如 `Chr1`、`So120_Chr1`），避免 `Chr_1`、`scaffold_1`。可用 chr_rename / fasta_id_renamer 批量改。

**Q3：`--step` 能只跑某一步吗？**
不能。CentIER 原版不支持单步执行，该参数仅作兼容保留，实际始终跑完整流程。

**Q4：会断点续传吗？**
Hi-C 自动模式的 HiC-Pro 步骤支持跳过已完成产物（除非 `--force-hicpro`）；CentIER 鉴定步骤本身每次重跑。想复用已有 Hi-C 矩阵可直接用 `--matrix1/2 --bed1/2` 传入。

**Q5：CentIER 返回非零退出码但生成了结果？**
程序会检查输出文件：只要关键产物（如 `centromere_range.txt`）生成就视为成功，仅在既失败又无产物时才报错。
