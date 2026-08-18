# 蛋白质到基因组比对 | Protein to Genome Alignment

一句话理解：**把一批蛋白质序列「对回」基因组上，找出每条蛋白落在基因组哪个位置、内含子/外显子结构如何，并顺手导出 GFF3/BED 注释和统计报告。**

## 功能概述 | Overview { #overview }

- 用 Miniprot 把蛋白质比对到基因组：擅长跨越内含子、处理远缘蛋白，是「蛋白到基因组」对齐的主流工具
- 自动解析 PAF 结果，统计序列一致性、覆盖度、mapping quality 等指标
- 一键导出多种产物：GFF3 注释(protein_match)、BED 坐标、按比对区间提取的基因组序列
- 单命令全流程：比对 → 解析 → 统计 → 导出注释 → 提取序列

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools pep2genome --genome genome.fa --protein protein.fa -o results
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 蛋白质 FASTA | 由 20 种氨基酸字母拼成的序列文件，本工具的「查询」 |
| 基因组 FASTA | 由 A/T/C/G 拼成的染色体序列文件，本工具的「地图」 |
| 比对(alignment) | 把两条序列「对齐」，看哪些位置相同、哪里断开 |
| 内含子 / 外显子 | 基因里最终被剪掉(内含子) / 保留(外显子)的部分 |
| 剪接(splicing) | 转录后把内含子去掉、外显子拼起来的过程 |
| PAF | 记录每条比对结果的制表符分隔文本格式 |
| 一致性(identity) | 两条序列相同位置的比例，越高越像 |
| 覆盖度(coverage) | 比对覆盖了整条序列的多大比例 |
| mapping quality | 比对置信度打分，越高越可靠 |
| GFF3 / BED | 两种记录「基因组上某区间」的标准坐标格式 |

## 输入 | Input { #input }

两个 FASTA 文件：

- `--genome`：基因组 FASTA，序列名任意
- `--protein`：蛋白质 FASTA，每条序列是一条查询蛋白

```text
>gene1 protein
MSKTRKVL...
>gene2 protein
MAAQ...
```

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required { #parameters-required }

**通俗理解|In plain words:** 三个必填：基因组「地图」、蛋白「查询」、输出目录。没有默认值，必须显式给出。

### 计算资源 | Compute { #parameters-compute }

**通俗理解|In plain words:** `--threads` 控制比对的并行度，越大越快（受机器 CPU 核数限制），一般不用动。`--miniprot-path` 是 miniprot 程序的位置，默认自动解析 annot 功能域环境，缺失时回退系统 PATH 查找，只有把 miniprot 装在非标准位置时才需要手动指定。

### 输出开关 | Output toggles { #parameters-output }

**通俗理解|In plain words:** 默认所有产物全出（PAF、统计、GFF3、BED、提取序列）。`--no-gff3` / `--no-bed` / `--no-statistics` 分别关掉对应产物，不需要哪样就关哪样，能省一点时间。

## 分析流程 | Pipeline { #pipeline }

```text
输入 genome.fa + protein.fa
    │
    ▼
miniprot 比对  ──►  alignment.paf
    │
    ▼
解析 PAF(一致性 / 覆盖度 / mapping quality)
    │
    ▼
统计报告  ──►  alignment_statistics.txt
    │
    ▼
导出 GFF3 / BED
    │
    ▼
seqkit subseq 按 BED 提取基因组序列  ──►  alignment.fa
```

## 输出 | Output { #output }

```text
results/
├── alignment.paf                     # miniprot 原始比对结果
├── alignment_statistics.txt          # 统计报告
├── alignment.gff3                    # GFF3 注释(protein_match)
├── alignment.bed                     # BED6 坐标
├── alignment.fa                      # 按 BED 提取的基因组序列
└── pep2genome_processing_YYYYMMDD_HHMMSS.log   # 运行日志
```

- `alignment.paf`：miniprot 的原始输出，一行一条比对，前 12 列是标准 PAF 字段
- `alignment_statistics.txt`：总体 / 查询 / 目标 / 质量 / 覆盖度五组统计
- `alignment.gff3`：每条蛋白比对成一个 `protein_match` 特征，属性里带 Identity / 覆盖率 / 比对得分
- `alignment.bed`：BED6(染色体、起、止、名字、分数、链)，可直接喂给 IGV 或 seqkit
- `alignment.fa`：seqkit 按 BED 从基因组切出的比对区序列

## 结果解读 | Interpreting Results { #interpreting-results }

- 统计报告：重点看「高质量比对数(MQ≥30 且 Identity≥80%)」占比，越高说明比对越可信
- GFF3 / BED：一条蛋白可能对应多条比对记录(同源多拷贝)，同一蛋白的多条记录可按 mapping quality 挑最佳
- PAF 为空：说明没有任何蛋白比对上，检查蛋白与基因组是否同源、序列 ID 是否正常
- 一致性(Identity)：反映蛋白与基因组翻译产物的相似程度，远缘物种会偏低

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 默认全开即可，最常见的需求就是「拿 GFF3 去注释基因组」
- 只需要坐标不要序列时，直接忽略 `alignment.fa` 即可
- 机器核数多可适当调大 `--threads` 加速 miniprot
- 基因组很大(如植物)时比对耗时较长，耐心等待(内部 24 小时超时)

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome, -g` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--protein, -p` | 必填 |  | 蛋白质FASTA文件｜Protein FASTA file |
| `--output-dir, -o` | `./pep2genome_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--miniprot-path` | `miniprot` |  | Miniprot工具路径｜Miniprot tool path |
| `--no-gff3` | — |  | 不导出GFF3格式｜Do not export GFF3 format |
| `--no-bed` | — |  | 不导出BED格式｜Do not export BED format |
| `--no-statistics` | — |  | 不生成统计报告｜Do not generate statistics report |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--protein` | 必填 |  | 蛋白质FASTA文件｜Protein FASTA file |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--miniprot-path` | — |  | Miniprot工具路径(默认域环境自动解析)｜Miniprot tool path (default: auto domain env) |
| `--no-gff3` | — | store_false | 不导出GFF3格式｜Do not export GFF3 format |
| `--no-bed` | — | store_false | 不导出BED格式｜Do not export BED format |
| `--no-statistics` | — | store_false | 不生成统计报告｜Do not generate statistics report |
| `--no-extract` | — | store_false | 不提取基因组序列｜Do not extract genome sequences |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- miniprot：蛋白质到基因组比对核心，自动解析 annot 域环境并经 conda run 调用，可用 --miniprot-path 或环境变量 MINIPROT_PATH 覆盖；域环境缺失时回退 PATH 直接调用
- seqkit：按 BED 提取基因组序列，自动解析 misc 域环境并经 conda run 调用，可用环境变量 SEQKIT_PATH 覆盖；域环境缺失时回退 PATH 直接调用
- Python 3 + numpy：统计计算

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
不支持。重新运行会覆盖同名输出文件(alignment.paf 等)。换参数重跑请换输出目录。

**Q2：报错 miniprot 找不到？**
用 `--miniprot-path /完整/路径/miniprot` 指定，或把 miniprot 加入 PATH。

**Q3：seqkit 缺失会怎样？**
比对、统计、GFF3 / BED 都正常，只有最后提取序列(alignment.fa)这一步失败。不需要序列时可直接忽略。

**Q4：为什么有的蛋白有多条 GFF3 记录？**
同源基因多拷贝时，一条蛋白会比对到基因组多个位置，每个位置对应一条记录，属正常现象。

**Q5：如何只看最佳比对？**
PAF 里每条记录有 mapping quality 和 identity，按 mapping quality 优先、identity 次之排序取第一条即可。