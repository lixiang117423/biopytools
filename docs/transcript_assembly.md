# 转录本从头组装分析模块 | Transcript De Novo Assembly Module

**基于参考基因组的转录本从头组装工具 | De Novo Transcript Assembly Tool Based on Reference Genome**

## 功能概述 | Overview

转录本从头组装分析模块是一个基于HISAT2比对和StringTie组装的转录本从头组装工具。适用于无GFF注释文件的场景，通过二代RNA-seq双端测序数据，对参考基因组进行转录本预测和组装，最终输出合并的转录本cDNA序列。支持多样本并行处理和断点续传。

## 主要特性 | Key Features

- **完整组装流程**: HISAT2索引 -> 比对(--dta) -> BAM排序 -> StringTie组装 -> GTF合并 -> 转录本提取
- **无需GFF注释**: 基于参考基因组直接从头组装转录本，不依赖已有注释
- **多样本支持**: 自动解析输入目录中的多样本FASTQ文件对
- **断点续传**: 已完成的步骤自动跳过，支持从中断处继续
- **灵活步骤控制**: 支持6个步骤的独立运行或完整流程
- **Conda环境支持**: 自动检测conda环境中的软件

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 完整的转录本从头组装流程
biopytools transcript-assembly \
    -g genome.fasta \
    -i ./clean_data \
    -o ./transcript_output
```

### 指定线程数和FASTQ模式 | Specify threads and FASTQ pattern

```bash
biopytools transcript-assembly \
    -g genome.fasta \
    -i ./clean_data \
    -o ./transcript_output \
    -t 24 \
    -p "*_1.clean.fq.gz"
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-g, --genome` | 参考基因组FASTA文件路径 | `-g genome.fasta` |
| `-i, --input` | 输入FASTQ文件目录 | `-i ./clean_data` |
| `-o, --output` | 输出目录 | `-o ./output` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-p, --pattern` | `*_1.fq.gz` | FASTQ文件命名模式（*为样本名占位符） |
| `-t, --threads` | `16` | 线程数 |
| `--sample-timeout` | `43200` | 单个样本处理超时时间（秒），默认12小时 |

### 步骤控制 | Step Control

| 参数 | 描述 |
|------|------|
| `-s, --step 1` | 构建HISAT2基因组索引 |
| `-s, --step 2` | HISAT2比对（开启--dta） |
| `-s, --step 3` | SAM转排序BAM |
| `-s, --step 4` | StringTie逐样本从头组装 |
| `-s, --step 5` | StringTie合并所有样本GTF |
| `-s, --step 6` | gffread提取转录本cDNA序列 |

### 高级选项 | Advanced Options

| 参数 | 描述 |
|------|------|
| `--force` | 强制重新处理已完成的步骤 |
| `-v` | 增加输出详细程度（可多次使用） |
| `--quiet` | 静默模式，仅输出错误信息 |

## 输入文件格式 | Input File Formats

### 参考基因组 | Reference Genome

标准FASTA格式的参考基因组序列（无需GFF注释文件）：

```
>chromosome1
ATCGATCGATCGATCGATCGATCG...
>chromosome2
GCTAGCTAGCTAGCTAGCTAGCTA...
```

### 测序数据 | Sequencing Data

存放于同一目录下的双端测序FASTQ文件，默认命名模式为 `*_1.fq.gz` / `*_2.fq.gz`：

```
clean_data/
├── sample1_1.fq.gz
├── sample1_2.fq.gz
├── sample2_1.fq.gz
├── sample2_2.fq.gz
└── sample3_1.fq.gz
└── sample3_2.fq.gz
```

## 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
transcript_output/
├── 01_hisat2_index/          # HISAT2基因组索引
│   └── genome.ht2*           # 索引文件
├── 02_hisat2_align/          # HISAT2比对结果
│   ├── sample1.sam
│   ├── sample2.sam
│   └── sample3.sam
├── 03_bam_sort/              # 排序后的BAM文件
│   ├── sample1.sorted.bam
│   ├── sample1.sorted.bam.bai
│   └── ...
├── 04_stringtie/             # 逐样本组装GTF
│   ├── sample1.gtf
│   ├── sample2.gtf
│   ├── sample3.gtf
│   └── gtf_list.txt          # GTF列表文件
├── 05_merge/                 # 合并后的GTF
│   └── merged.gtf
├── 06_transcripts/           # 最终转录本序列
│   └── transcripts.fa        # 转录本cDNA序列
└── 99_logs/                  # 日志文件
    └── pipeline.log
```

### 关键输出文件 | Key Output Files

| 文件 | 描述 |
|------|------|
| `05_merge/merged.gtf` | 所有样本合并后的统一GTF注释 |
| `06_transcripts/transcripts.fa` | 最终输出的转录本cDNA序列（FASTA格式） |

## 流程说明 | Pipeline Details

### Step 1: HISAT2索引构建

使用`hisat2-build`对参考基因组构建索引。索引文件存储在`01_hisat2_index/`目录下。

### Step 2: HISAT2比对

使用`hisat2`将每个样本的clean reads比对到参考基因组。开启`--dta`参数以优化后续StringTie的转录本组装效果。

### Step 3: SAM转BAM

使用`samtools sort`将比对结果从SAM格式转换为排序后的BAM格式，并使用`samtools index`建立索引。

### Step 4: StringTie逐样本组装

使用`stringtie`对每个样本独立进行转录本从头组装。不提供`-G`参数（参考GTF），实现真正的从头转录本预测。

### Step 5: StringTie合并

使用`stringtie --merge`将所有样本的GTF文件合并为一个统一的`merged.gtf`。合并时需要基因组FASTA作为参考。

### Step 6: 转录本序列提取

使用`gffread`根据`merged.gtf`和基因组FASTA文件，提取所有转录本的cDNA序列，输出为`transcripts.fa`。

## 依赖软件 | Dependencies

- **HISAT2** - 基因组比对
- **SAMtools** - BAM文件处理
- **StringTie** - 转录本组装和合并
- **gffread** (gffutils) - 转录本序列提取

建议使用Conda环境安装上述软件。

## 常见问题 | Troubleshooting

**Q: 找不到样本FASTQ文件对**

检查输入目录中的文件命名是否符合默认模式`*_1.fq.gz`/`*_2.fq.gz`。如果不是，请使用`-p`参数指定正确的命名模式，例如：`-p "*_R1.fq.gz"`。

**Q: HISAT2比对失败**

确认参考基因组FASTA文件格式正确。检查线程数是否超出系统限制。

**Q: StringTie合并步骤失败**

确保所有样本都成功完成了Step 4的组装。检查`04_stringtie/`目录下是否存在所有样本的GTF文件。

**Q: 如何从中断处继续运行**

直接重新运行相同的命令即可。已完成的步骤会自动跳过（断点续传）。如需强制重新运行所有步骤，请添加`--force`参数。
