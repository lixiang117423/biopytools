# insert-detection - 插入序列位点检测 | Insert Sequence Insertion Site Detection

一句话理解：**从配对测序数据（FASTQ）里找出某段外源插入序列（如 T-DNA）插到了参考基因组的哪个位置**——先把 reads 比对到基因组，再把「一半比对在基因组、一半悬空（soft-clip）」的 reads 的悬空部分比对回插入序列，从而锁定插入位点并打分。

## 功能概述 | Overview

- 检测插入序列（T-DNA 等）在参考基因组上的插入位点，输出每个位点的坐标、方向、支持 reads 数与得分
- 三步流程：reads 比对基因组 → 提取 soft-clip 序列 → soft-clip 序列回比插入序列，再按 5 组 reads 分类加权打分
- 自动识别 FASTQ 目录中的成对样品（默认匹配 biopytools fastp 输出后缀 _1.clean.fq.gz / _2.clean.fq.gz，可自定义）
- 断点续传：每步按输出文件存在性跳过已完成部分，--force 强制全部重跑
- 自动构建 bowtie2 索引（基因组和插入序列各一次），生成可导入 IGV 的 BAM 方便人工核对

## 快速开始 | Quick Start

```bash
biopytools insert-detection -i genome.fa --insert tdna.fa --fastq-dir fastq_output/ -o output/
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 参考基因组 genome | 物种本来的 DNA 序列（FASTA），作为「地图」 |
| 插入序列 insert | 外源插入的那段序列（如 T-DNA），作为「要找的标记」 |
| 配对测序 FASTQ | 测序仪读出的短片段，每个样品分 R1/R2 两个文件 |
| soft-clip 软裁剪 | 一条 read 比对时只有一部分能对上基因组，悬空的另一部分被「裁剪」出来，正是跨越插入位点边界的那半截 |
| CIGAR | 比对结果里描述 read 哪段比对上、哪段被裁剪的字符串（S 表示 soft-clip） |
| 支持 reads | 指向同一个插入位点的 read 数量，越多越可信 |
| 得分 score | 综合 5 类 read 证据加权算出的可信度，超过阈值才算候选位点 |
| BAM / SAM | 比对结果的标准二进制/文本格式 |

## 输入 | Input

### 参考基因组 FASTA

-i/--genome 指定，FASTA 格式。程序会自动构建 bowtie2 索引（若 .fa 去掉后缀后的索引文件不存在）。

### 插入序列 FASTA

--insert 指定，FASTA 格式（支持 .fa/.fasta/.fna），程序自动构建索引。

### FASTQ 目录

--fastq-dir 指定一个目录，程序扫描其中成对的文件。默认匹配 fastp 输出命名：R1 后缀 _1.clean.fq.gz、R2 后缀 _2.clean.fq.gz，可用 --read1-suffix/--read2-suffix 自定义。

```text
fastq_output/
├── sampleA_1.clean.fq.gz
├── sampleA_2.clean.fq.gz
├── sampleB_1.clean.fq.gz
└── sampleB_2.clean.fq.gz
```

## 参数说明 | Parameters

### 输入与输出 | Input & output

**通俗理解|In plain words:** -i 基因组、--insert 插入序列、--fastq-dir FASTQ 目录三个必填，-o 输出目录（默认 ./insert_detection_output）。

### 检测参数 | Detection

**通俗理解|In plain words:** --min-clip 是最小 soft-clip 长度（低于此长度的悬空段不算证据，调大更严格）；--score-threshold 是得分阈值（调大只留更可信的位点、可能漏掉弱信号，调小更敏感但假阳性增多）。注意 --min-support（最小支持 reads 数）目前只是声明并记录，实际过滤靠得分阈值。一般用默认即可。

### 工具路径 | Tool paths

**通俗理解|In plain words:** --bowtie2-path、--samtools-path 指定软件路径，默认自动解析功能域环境，缺失时回退 PATH 里的 bowtie2、samtools。一般不用动，除非软件不在 PATH。

### FASTQ 命名与流程控制 | File pattern & control

**通俗理解|In plain words:** --read1-suffix/--read2-suffix 告诉程序「R1/R2 文件名长什么样」，只有文件后缀匹配才会被当成样品；--force 强制全部重跑（覆盖断点续传）；--verbose/--quiet 控制日志详细程度。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 核心思想是「跨界 read 的两头各比对一次」——先比对基因组看它落在哪，再把悬空的 soft-clip 部分比对回插入序列，从而同时确定基因组位点和插入序列位点。

```text
输入: 基因组 FASTA + 插入序列 FASTA + FASTQ 目录
    |
    v
识别成对样品 (按 R1/R2 后缀匹配)
    |
    v
步骤1: reads 比对到基因组 (bowtie2 --very-sensitive-local, 生成 sorted.bam)
    |
    v
步骤2: 提取 soft-clip 序列 (samtools + 自定义脚本, 生成 softclip.fastq + 供 IGV 的 bam)
    |
    v
步骤3: soft-clip 序列比对到插入序列 (bowtie2 --local -L 10, 生成 tdna.bam)
    |
    v
步骤4: 过滤与评分 (错配<=1, 按边界远近设匹配长度阈值, 5 组 reads 加权打分)
    |   (权重 700/300/300/1/1, 得分 >= score-threshold 才算候选位点)
    v
输出: insertion_sites.tsv + summary.txt
```

## 输出 | Output

```text
output/
├── 00_pipeline_info/
│   ├── software_versions.yml        # 软件版本
│   └── pipeline_params.yml          # 运行参数
├── 04_results/
│   ├── insertion_sites.tsv         # 插入位点总表 (核心)
│   └── summary.txt                  # 汇总统计
├── 99_logs/
│   └── insert_detection.log         # 运行日志
├── tmp/                             # 临时文件 (运行结束清理)
└── <sample>_insert_detection/       # 每个样品一个目录
    ├── 01_genome_alignment/<sample>.sorted.bam (+.bai)
    ├── 02_softclip_extraction/<sample>.softclip.fastq (+.softclip.bam)
    └── 03_insert_alignment/<sample>.tdna.bam (+.bai)
```

### insertion_sites.tsv 字段

| 字段 | 说明 |
|------|------|
| sample_id | 样品名 |
| chromosome / start / end | 基因组上的插入位点坐标 |
| orientation | 插入方向（+/-） |
| insert_position | 插入序列上的断点位置（tdna::pos） |
| junction | 连接类型（T-DNA:REF 或 REF:T-DNA） |
| support_reads | 支持该位点的 reads 总数 |
| score | 加权得分（越高越可信） |
| group_counts | 5 组 reads 各自计数（逗号分隔） |

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** insertion_sites.tsv 是最终答案——每个候选插入位点一行，得分越高、支持 reads 越多，越可信。

- 优先看 score 高的位点：得分由 5 组 reads 加权（700/300/300/1/1），分值越高证据越强
- support_reads 越多越可信；单条 reads 支撑的位点（group_counts 里大部分组为 0）要谨慎
- orientation 与 junction 描述插入方向：junction 为 T-DNA:REF 表示断点在 T-DNA 左侧，REF:T-DNA 表示右侧
- summary.txt 给出每个样品的位点数、染色体分布、得分范围，可快速判断哪些样品有插入事件
- 用 03_insert_alignment 下的 tdna.bam 在 IGV 里人工复核高分位点的比对情况

## 参数选择建议 | Parameter Guidance

- 常规检测：全部用默认即可，重点看 insertion_sites.tsv 的 score 与 support_reads
- 假阳性偏多：调高 --score-threshold（如 2000）或调大 --min-clip（如 30），只留更可信的位点
- 弱信号漏检：调低 --score-threshold 或调小 --min-clip，但假阳性会增多
- FASTQ 不是 fastp 输出命名：用 --read1-suffix/--read2-suffix 指定实际后缀
- 换参数重跑：加 --force，否则断点续传会复用旧结果

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `--insert` | 必填 |  | 插入序列FASTA文件｜Insert sequence FASTA file |
| `--fastq-dir` | 必填 |  | FASTQ文件目录｜FASTQ files directory |
| `-o, --output-dir` | `./insert_detection_output` |  | 输出目录｜Output directory (default: ./insert_detection_output) |
| `-t, --threads` | `12` |  | 线程数｜Threads |
| `--min-clip` | `20` | int | 最小soft-clip长度｜Minimum soft-clip length |
| `--min-support` | `5` | int | 最小支持reads数｜Minimum supporting reads |
| `--score-threshold` | `1000` | int | 得分阈值｜Score threshold |
| `--bowtie2-path` | `bowtie2` |  | Bowtie2可执行文件路径｜Bowtie2 executable path |
| `--samtools-path` | `samtools` |  | samtools可执行文件路径｜samtools executable path |
| `--read1-suffix` | `_1.clean.fq.gz` |  | R1文件后缀（包含扩展名）｜Read 1 file suffix with extension |
| `--read2-suffix` | `_2.clean.fq.gz` |  | R2文件后缀（包含扩展名）｜Read 2 file suffix with extension |
| `--force` | — |  | 强制重新运行所有步骤｜Force rerun all steps |
| `--verbose` | — |  | 显示详细日志｜Verbose logging |
| `--quiet` | — |  | 仅显示错误｜Errors only |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `--insert` | 必填 |  | 插入序列FASTA文件｜Insert sequence FASTA file |
| `--fastq-dir` | 必填 |  | FASTQ文件目录｜FASTQ files directory |
| `-o, --output-dir` | `./insert_detection_output` |  | 输出目录｜Output directory (default: ./insert_detection_output) |
| `-t, --threads` | `12` | int | 线程数｜Threads (default: 12) |
| `--min-clip` | `20` | int | 最小soft-clip长度｜Minimum soft-clip length (default: 20) |
| `--min-support` | `5` | int | 最小支持reads数｜Minimum supporting reads (default: 5) |
| `--score-threshold` | `1000` | int | 得分阈值｜Score threshold (default: 1000) |
| `--bowtie2-path` | — |  | Bowtie2可执行文件路径(默认域环境自动解析)｜Bowtie2 executable path (default: auto domain env) |
| `--samtools-path` | — |  | samtools可执行文件路径(默认域环境自动解析)｜samtools executable path (default: auto domain env) |
| `--read1-suffix` | `_1.clean.fq.gz` |  | R1文件后缀（包含扩展名，默认匹配fastp输出）｜Read 1 file suffix with extension (default: _1.clean.fq.gz, matches fastp output) |
| `--read2-suffix` | `_2.clean.fq.gz` |  | R2文件后缀（包含扩展名，默认匹配fastp输出）｜Read 2 file suffix with extension (default: _2.clean.fq.gz, matches fastp output) |
| `--force` | — | store_true | 强制重新运行所有步骤｜Force rerun all steps |
| `--verbose` | — | store_true | 显示详细日志｜Show verbose logs |
| `--quiet` | — | store_true | 仅显示错误日志｜Show error logs only |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3（需 pandas、numpy）
- bowtie2（bowtie2、bowtie2-build，无对应功能域环境，默认从 PATH 调用，可用 --bowtie2-path 或环境变量 BOWTIE2_PATH / BOWTIE2_BUILD_PATH 覆盖）
- samtools（自动解析 align 域环境并经 conda run 调用，可用 --samtools-path 或环境变量 SAMTOOLS_PATH 覆盖；域环境缺失时回退 PATH 直接调用）
- minimap2（可选，用于长序列；自动解析 align 域环境，可用环境变量 MINIMAP2_PATH 覆盖）

## 常见问题 | FAQ

**Q1：能断点续传吗？**
能。基因组比对、soft-clip 提取、插入序列比对各步都按输出文件存在性跳过已完成部分；用 --force 强制全部重跑。

**Q2：--min-support 设置了为什么好像没效果？**
--min-support（最小支持 reads 数）目前只被声明和记录到 pipeline_params.yml，实际过滤靠 --score-threshold 的加权得分。若想收紧，请调 --score-threshold 而不是 --min-support。

**Q3：为什么没识别到我的 FASTQ 样品？**
程序只认文件名后缀匹配 --read1-suffix/--read2-suffix 的文件（默认 _1.clean.fq.gz/_2.clean.fq.gz），且必须 R1/R2 成对出现。请确认命名或自定义后缀。

**Q4：索引会自动构建吗？**
会。基因组和插入序列的 bowtie2 索引若不存在，程序会自动运行 bowtie2-build 构建。

**Q5：soft-clip 里靠近边界和远离边界的阈值不一样吗？**
是。比对到插入序列时，靠近边界的 reads 匹配长度达到 18 即可，远离边界的需达到 30（内部固定，未暴露为参数），这是为了平衡边界证据与内部证据的可靠性。