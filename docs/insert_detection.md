# insert-detection | 插入序列位点检测

**基于 Bowtie2 + soft-clip read 信号，从重测序数据中定位外源插入序列（如 T-DNA、转座子）在宿主基因组中的插入位点 | Detect insertion sites of an exogenous insert sequence in a host genome from paired-end resequencing reads**

## 功能概述 | Overview

本模块面向功能基因组学中常见的"插入序列定位"问题：已知一段外源序列（例如 T-DNA 载体、转座子标签、报告基因），需要从携带该插入的个体重测序数据中找到它在参考基因组上的精确插入位置。

核心策略是：用 Bowtie2 将 read 分别比对到"参考基因组 + 插入序列"组合索引上，从比对结果中寻找一边唯一比对到插入序列、另一端 soft-clip 到参考基因组的嵌合 read（split-read），根据这些 read 的剪接点对候选插入位点进行打分聚类，最终输出候选插入位点、支持 read 数与得分。流程自动识别 FASTQ 目录中的配对样本（默认匹配 fastp 输出的 `_1.clean.fq.gz` / `_2.clean.fq.gz`）。

## 快速开始 | Quick Start

```bash
biopytools insert-detection \
    -i host_genome.fa \
    --insert tdna_vector.fa \
    --fastq-dir ./fastq_output/ \
    -o ./insert_detection_output/
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --genome` | 参考基因组 FASTA 文件 |
| `--insert` | 插入序列 FASTA 文件 |
| `--fastq-dir` | FASTQ 文件目录（自动识别配对样本） |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./insert_detection_output` | 输出目录 |
| `-t, --threads` | `12` | 线程数 |
| `--min-clip` | `20` | 最小 soft-clip 长度 |
| `--min-support` | `5` | 最小支持 read 数 |
| `--score-threshold` | `1000` | 得分阈值 |
| `--bowtie2-path` | `bowtie2` | Bowtie2 可执行文件路径 |
| `--samtools-path` | `samtools` | samtools 可执行文件路径 |
| `--read1-suffix` | `_1.clean.fq.gz` | R1 文件后缀（含扩展名） |
| `--read2-suffix` | `_2.clean.fq.gz` | R2 文件后缀（含扩展名） |
| `--force` | 关闭 | 强制重新运行所有步骤 |
| `--verbose` | 关闭 | 显示详细日志 |
| `--quiet` | 关闭 | 仅显示错误 |

（运行 `biopytools insert-detection -h` 查看完整参数列表）

## 输出 | Output

```
insert_detection_output/
├── 04_results/
│   ├── insert_sites.tsv          # 候选插入位点（染色体、位置、链、支持 read 数、得分）
│   └── insert_summary.txt        # 运行总结
└── 99_logs/                      # 每个样本/步骤的日志
```

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
| `--bowtie2-path` | `bowtie2` |  | Bowtie2可执行文件路径｜Bowtie2 executable path |
| `--samtools-path` | `samtools` |  | samtools可执行文件路径｜samtools executable path |
| `--read1-suffix` | `_1.clean.fq.gz` |  | R1文件后缀（包含扩展名，默认匹配fastp输出）｜Read 1 file suffix with extension (default: _1.clean.fq.gz, matches fastp output) |
| `--read2-suffix` | `_2.clean.fq.gz` |  | R2文件后缀（包含扩展名，默认匹配fastp输出）｜Read 2 file suffix with extension (default: _2.clean.fq.gz, matches fastp output) |
| `--force` | — | store_true | 强制重新运行所有步骤｜Force rerun all steps |
| `--verbose` | — | store_true | 显示详细日志｜Show verbose logs |
| `--quiet` | — | store_true | 仅显示错误日志｜Show error logs only |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)（必需，比对）
- [samtools](http://www.htslib.org/)（必需，BAM 排序与索引）
- Python：标准库

## 引用 | Citation

- Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. *Nature Methods*, 2012, 9(4): 357-359. doi:10.1038/nmeth.1923
- Li H et al. The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 2009, 25(16): 2078-2079. doi:10.1093/bioinformatics/btp352

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [Bowtie2 手册](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [samtools 文档](http://www.htslib.org/doc/samtools.html)
