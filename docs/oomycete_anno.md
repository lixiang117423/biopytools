# oomycete_anno — 疫霉菌基因组注释 | Oomycete Genome Annotation (T2T Augustus)

一句话理解：**给疫霉菌（卵菌）基因组自动「写基因说明书」**——复刻 T2T 证据驱动的 Augustus 流程，用 RNA-seq / 同源蛋白 / 三代转录本等证据辅助，从头预测基因结构，输出一份 GFF3 注释，并对已知效应子位点做专门「救援」。

## 功能概述 | Overview

- 复刻 T2T 证据驱动 Augustus 手工流程（Phase1 + Phase2 全量）
- 九步流水线：重复屏蔽 → RNA-seq 比对 → 三代转录本 / 蛋白比对 → hints → 训练集 → Augustus 训练 + 预测 → 效应子救援 → LTR 注释
- 证据 graceful degradation：缺某类证据自动跳过对应步骤，仍能 ab initio 预测
- Phase3 效应子救援：用已知效应子的全长 miniprot 比对当基因模型，替换/补回 Augustus 在效应子簇的错注漏注
- 断点续传：每步按关键输出文件存在性跳过已完成部分
- 自动处理中文路径坑（C 二进制与 GeneMark 的中文路径乱码），把计算放到 `~/tmp` 的 ASCII 临时目录

## 快速开始 | Quick Start

```bash
biopytools oomycete-anno -g genome.fa -s phytophthora --rnaseq-dirs rna1/ rna2/ -o out/ -t 24
```

最小输入：一个基因组 FASTA + 一个 Augustus 物种名。RNA-seq 等证据都是可选的，没有也能跑（纯 ab initio）。

## 零基础概念速览 | Concepts in plain words

不熟悉生信术语的话，先花两分钟看这张表：

| 术语 | 通俗理解 |
|------|----------|
| 基因注释 | 在基因组上标出「哪里有基因、外显子怎么切、编码什么」，像给一长串 DNA 写目录 |
| 从头预测 (ab initio) | 不靠已知基因、只凭序列本身的规律猜基因结构，像「盲听写」 |
| 证据 (evidence) | RNA-seq 等实验数据，用来「提示」基因在哪、边界在哪，让预测更准 |
| Augustus | 一个经典基因预测软件，会用训练集学到的规律 + hints 证据来预测 |
| 训练集 (training set) | 用来「教」Augustus 认识这个物种基因特征的一组已知基因模型 |
| hints | 从 RNA-seq/蛋白比对提取的「线索」（如剪接位点、编码区），作为软证据喂给 Augustus |
| 重复序列 / 重复屏蔽 | 基因组里大量重复的「噪音」片段；先软屏蔽（打成小写）避免干扰基因预测 |
| 效应子 (effector) | 病原菌分泌出来攻击宿主的小蛋白，卵菌致病的关键武器 |
| T2T 流程 | 「端到端」的高质量注释流程，本文特指证据驱动的 Augustus 注释法 |

## 输入 | Input

### 基因组 FASTA（-g, --genome）

待注释的疫霉菌基因组，标准 FASTA。

### Augustus 物种名（-s, --species）

给 Augustus 训练出来的物种参数命名，用简单字母（如 `phytophthora`）。

### 可选证据 | Optional evidence

- `--rnaseq-dirs`：二代 RNA-seq 目录（可多个），文件按 R1/R2 后缀配对（见下）
- `--prot-seq`：同源蛋白文件（Phase2 蛋白 hints）
- `--isoseq`：三代转录本文件（Phase2，供 TransDecoder 生成训练集）
- `--effectors`：已知效应子蛋白（Phase3 救援）

RNA-seq 配对默认按 `_1.clean.fq.gz` / `_2.clean.fq.gz` 识别，可用 `--read1-pattern` / `--read2-pattern` 调整。

## 参数说明 | Parameters

### 必需与证据 | Required & evidence

**通俗理解|In plain words:** `-g` 和 `-s` 必填。四类证据（RNA-seq、同源蛋白、三代转录本、效应子）都是「可选项」——给得越多注释越准，但全不给也能跑纯 ab initio。链特异性 RNA-seq 数据要用 `--rna-strandness` 标注 FR/RF，普通非链特异性数据留空即可。

### 流程开关 | Pipeline skip flags

**通俗理解|In plain words:** 这一组 `--skip-*` 开关用于「跳过某个证据步骤」。比如没有三代转录本就别传 `--isoseq`（自动跳过），或者已手动做完了重复屏蔽想复用就用 `--skip-repeat`。**正常情况下一个都不用动，程序会根据你有没有提供对应证据自动判断。**

### 运行参数 | Runtime

**通俗理解|In plain words:** `-t` 是线程数；`--no-soft-masking` 把重复屏蔽从软屏蔽（小写）改成硬屏蔽（N），一般保持默认软屏蔽；`--gmes-petap-path` / `--genemark-perl-env` 是 GeneMark 的安装路径与 perl 环境，装在不同位置才需要指定。

### 效应子救援参数 | Effector rescue

**通俗理解|In plain words:** `--rescue-min-identity` 是判定效应子「全长比对」的最低 identity（注意全长主要靠结构完整性判断，不是死卡高分）；`--rescue-conflict-overlap` 是「Augustus 基因与效应子模型重叠超过多少比例就算冲突、该替换」。**默认值经实践验证，一般不用动。**

## 分析流程 | Pipeline

**通俗理解|In plain words:** 九步流水线，缺哪类证据就跳过哪步，最后一定产出 Augustus 注释 GFF（有效应子时再救援一遍）。

```text
输入 基因组 + （可选）各类证据
    │
    ▼
01 重复屏蔽（RepeatModeler + RepeatMasker，软屏蔽）
    │
    ▼
02 RNA-seq 比对（HISAT2 -> sorted BAM）
    │
    ▼
03 三代转录本比对（GMAP -> 基因组坐标 GFF3，可选）
    │
    ▼
04 蛋白比对（miniprot -> CDSpart hints，可选）
    │
    ▼
05 hints（bam2hints 抽剪接位点 + 合并蛋白 hints）
    │
    ▼
06 训练集（有三代 -> TransDecoder；否则 GeneMark-ES）
    │
    ▼
07 Augustus 训练（etraining）+ 预测（augustus）
    │
    ▼
09 效应子救援（miniprot 全长比对替换/补回，可选，失败不阻断）
    │
    ▼
08 LTR 注释（gt ltrharvest + LTR_retriever，正交 TE，可选，失败不阻断）
```

## 输出 | Output

```text
out/
├── 00_pipeline_info/
│   └── software_versions.yml          # 软件版本与参数
├── 01_repeat_masking/                 # 屏蔽后的基因组 + RepeatModeler 库
├── 02_rna_align/                      # rnaseq.sorted.bam + 索引
├── 03_iso_align/                      # iso.gff3（三代转录本，可选）
├── 04_protein_align/                  # miniprot.gff3（蛋白，可选）
├── 05_hints/                          # intron.gff / protein.gff / hintsfile.gff
├── 06_training/                       # genemark.gtf 或 iso_training.gff3
├── 07_augustus/
│   └── {species}.augustus.gff         # 最终基因注释 GFF（核心）
├── 08_ltr/                            # LTR 注释产物（可选）
├── 09_effector_rescue/
│   └── {species}.rescued.gff          # 救援后的最终 GFF（核心，有 --effectors 时）
└── 99_logs/                           # 运行日志
```

### 关键文件说明 | Key files

- `{species}.augustus.gff`：Augustus 预测的基因注释，标准 GFF 格式，是主结果
- `{species}.rescued.gff`：提供 `--effectors` 且未 `--skip-rescue` 时生成，在 Augustus 结果基础上替换/补回了效应子位点的基因模型，**用它作最终注释**
- `hintsfile.gff`：合并后的证据 hints，可观察哪些位点有 RNA-seq / 蛋白证据支撑
- `software_versions.yml`：记录所有软件版本与运行参数，便于论文 Methods 复现

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 主结果是一份 GFF 注释文件，最终版本（有效应子时）是 rescue 后的那个。判断质量主要看基因数量是否合理、效应子位点是否被正确处理。

- **基因数量**：疫霉菌基因组通常编码 1–2 万左右的基因，数量明显偏少（几千）可能训练集不足或 hints 缺失，偏多（数万）可能是重复屏蔽不充分导致假基因
- **效应子位点**：卵菌效应子常成簇分布且彼此相似，Augustus 在此容易把相邻旁系同源基因合并成「嵌合基因」。救援步骤用 miniprot 逐条独立比对来纠正——若提供了 `--effectors`，优先检查 `09_effector_rescue` 里 rescue 替换/补回了多少基因（日志会打印「替换 N 个、加入 M 个」）
- **hints 覆盖**：`05_hints` 里 hints 越多，说明证据越充分，预测越可信；完全没有 hints 时 Augustus 是纯 ab initio，边界精度会下降

## 参数选择建议 | Parameter Guidance

- **有 RNA-seq**：给 `--rnaseq-dirs`，这是性价比最高的证据，能显著提升剪接边界精度
- **有高质量三代转录本**：给 `--isoseq`，会走 TransDecoder 生成训练集（T2T 正宗做法），比 GeneMark 训练更接近真实基因
- **做效应子研究**：务必给 `--effectors`，否则 Augustus 在效应子簇位点大概率错注
- **只想快速出个草稿**：只给 `-g -s`，纯 ab initio 跑（GeneMark-ES 自训练）
- **复用已有屏蔽**：手动做完了重复屏蔽就 `--skip-repeat`
- **链特异性数据**：据实验建库方向设 `--rna-strandness FR` 或 `RF`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组 FASTA｜Genome FASTA |
| `-s, --species` | 必填 |  | Augustus 物种名(简单字母, 如 phytophthora)｜Augustus species name |
| `--rnaseq-dirs` | — |  | 二代 RNA-seq 目录(可多个)｜Short RNA-seq dir(s) |
| `--prot-seq` | — |  | 同源蛋白文件(Phase2)｜Homologous proteins (P2) |
| `--isoseq` | — |  | 三代转录本文件(Phase2)｜Long-read transcripts (P2) |
| `--effectors` | — |  | 已知效应子蛋白(Phase3 救援)｜Known effectors (P3 rescue) |
| `--read1-pattern` | `_1.clean.fq.gz` |  | R1 文件后缀模式｜R1 suffix pattern |
| `--read2-pattern` | `_2.clean.fq.gz` |  | R2 文件后缀模式｜R2 suffix pattern |
| `--rna-strandness` | `` |  | 链特异性: ''(非链特异性)/FR/RF｜Strandness: ''/FR/RF |
| `-o, --output-dir` | `./oomycete_anno_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` |  | 线程数｜Threads |
| `--no-soft-masking` | `False` |  | 禁用软屏蔽(改用硬屏蔽)｜Disable soft masking |
| `--gmes-petap-path` | — |  | GeneMark gmes_petap.pl 路径｜GeneMark gmes_petap.pl path |
| `--genemark-perl-env` | — |  | GeneMark perl 提供环境｜GeneMark perl provider env |
| `--skip-repeat` | `False` |  | 跳过重复屏蔽｜Skip repeat masking |
| `--skip-rna` | `False` |  | 跳过 RNA-seq 比对｜Skip RNA-seq alignment |
| `--skip-iso` | `False` |  | 跳过三代转录本(Phase2)｜Skip long-read (P2) |
| `--skip-protein` | `False` |  | 跳过蛋白证据(Phase2)｜Skip protein hints (P2) |
| `--skip-ltr` | `False` |  | 跳过 LTR 注解(Phase2)｜Skip LTR annotation (P2) |
| `--skip-rescue` | `False` |  | 跳过效应子救援(Phase3)｜Skip effector rescue (P3) |
| `--rescue-min-identity` | `0.85` | float | 效应子救援 miniprot 最低 identity｜Rescue min identity |
| `--rescue-conflict-overlap` | `0.5` | float | Augustus 与效应子模型重叠>此比例则替换｜Conflict overlap fraction |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组 FASTA｜Genome FASTA |
| `-s, --species` | 必填 |  | Augustus 物种名(简单字母, 如 phytophthora)｜Augustus species name (simple, e.g. phytophthora) |
| `--rnaseq-dirs` | — |  | 二代 RNA-seq 目录(可多个)｜Short RNA-seq dir(s) |
| `--prot-seq` | — |  | 同源蛋白文件(Phase2)｜Homologous proteins (P2) |
| `--isoseq` | — |  | 三代转录本文件(Phase2)｜Long-read transcripts (P2) |
| `--effectors` | — |  | 已知效应子蛋白(Phase3 救援, 直接当基因模型替换错注/漏注位点)｜Known effectors (P3 rescue, used as gene models to fix loci) |
| `--read1-pattern` | `_1.clean.fq.gz` |  | R1 文件后缀模式｜R1 suffix pattern |
| `--read2-pattern` | `_2.clean.fq.gz` |  | R2 文件后缀模式｜R2 suffix pattern |
| `--rna-strandness` | `` |  | 链特异性: ''(非链特异性,默认) / FR / RF｜Strandness: ''(unstranded,default)/FR/RF |
| `-o, --output-dir` | `./oomycete_anno_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Threads |
| `--no-soft-masking` | — | store_false | 禁用软屏蔽(改用硬屏蔽)｜Disable soft masking (use hard masking) |
| `--gmes-petap-path` | — |  | GeneMark gmes_petap.pl 路径(默认 ~/software/GeneMark/...)｜GeneMark gmes_petap.pl path |
| `--genemark-perl-env` | — |  | GeneMark perl 提供环境(默认 braker_v.3.0.8)｜GeneMark perl provider env |
| `--skip-repeat` | — | store_true | 跳过重复屏蔽｜Skip repeat masking |
| `--skip-rna` | — | store_true | 跳过 RNA-seq 比对｜Skip RNA-seq alignment |
| `--skip-iso` | — | store_true | 跳过三代转录本(Phase2)｜Skip long-read (P2) |
| `--skip-protein` | — | store_true | 跳过蛋白证据(Phase2)｜Skip protein hints (P2) |
| `--skip-ltr` | — | store_true | 跳过 LTR 注解(Phase2)｜Skip LTR annotation (P2) |
| `--skip-rescue` | — | store_true | 跳过效应子救援(Phase3)｜Skip effector rescue (P3) |
| `--rescue-min-identity` | `0.85` | float | 效应子救援 miniprot 最低 identity(全长靠 Target起=1+stop 判)｜Rescue min identity |
| `--rescue-conflict-overlap` | `0.5` | float | Augustus 基因与效应子模型重叠>此比例则替换｜Conflict overlap fraction |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

以下软件通过 conda 环境自动检测调用（括号内为默认环境名）：

- RepeatModeler / RepeatMasker / BuildDatabase（`repeat` 环境）
- HISAT2 / hisat2-build（`rna` 环境）
- GMAP / gmap_build（`pasa_v.2.5.3` 环境）
- TransDecoder（`annot` 环境 + `transdecoder_v.5.5.0`）
- miniprot、bam2hints、augustus、etraining、gtf2gff、gff2gb、gt（`annot` / `Augustus_v.3.5.0` 环境）
- LTR_retriever（`repeat` 环境）
- GeneMark（`~/software/GeneMark/gmes_linux_64_4/gmes_petap.pl`，perl 环境 `braker_v.3.0.8`）

## 常见问题 | FAQ

**Q1：路径含中文，跑到 GeneMark / Augustus 就报乱码或失败？**
这是已知坑：GeneMark probuild、Augustus 等 C 二进制以及 gmes_petap.pl 会把中文路径解成乱码。程序已自动把计算放到 `~/tmp` 下的 ASCII 临时目录、运行完拷回产物，无需手动处理。

**Q2：为什么 GeneMark 只跑 ES 不跑 ET/EP+？**
非链特异性 RNA-seq 的 bam2hints 输出（score=0、strand=.）会让 GeneMark-ET 的 parse_ET.pl 除零崩溃（已实测）。因此程序始终用最稳健的 GeneMark-ES（基因组自训练），intron 与蛋白证据改走 Augustus hints，证据不浪费。

**Q3：bam2hints 报缺 libbamtools.so.2.5.2？**
Augustus 环境缺这个库，程序会自动从其它环境软链一个 ABI 兼容的 2.5.x 版到 `~/.cache/biopytools/lib` 兜底。

**Q4：为什么 Augustus 必须直接调二进制、不用 conda run？**
conda run 的激活会把 `AUGUSTUS_CONFIG_PATH` 覆盖成环境自带的 config，破坏程序准备好的物种 config，所以 Augustus/etraining 直接调二进制并手动注入环境变量。

**Q5：换参数重跑，为什么有些步骤没重新跑？**
断点续传按关键输出文件存在性跳过。要强制重跑某步，先删除对应步骤目录里的关键产物（如 `07_augustus/{species}.augustus.gff`），或整体换输出目录。

**Q6：效应子救援失败了，整个注释会失败吗？**
不会。救援和 LTR 都是「失败不阻断」的可选步骤，失败会打警告并保留 Augustus 结果。
