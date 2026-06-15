# Bismark甲基化分析 | Bismark Methylation Analysis

**基于Bismark的全基因组重亚硫酸盐测序（WGBS）甲基化分析流程 | Whole-genome bisulfite sequencing methylation analysis pipeline based on Bismark**

## 功能概述 | Overview

封装 Bismark + Bowtie2 工具链，对一批已质控的配对端 FASTQ 进行重亚硫酸盐比对、甲基化提取与汇总。流程会自动构建基因组索引、按 R1 文件后缀模式自动匹配 R2、调用 `bismark` 比对与 `bismark_methylation_extractor` 提取 CPG/CHG/CHH 甲基化信息，并生成最终汇总报告。

- 自动完成 Bismark 基因组索引构建（已存在则跳过）
- 按文件后缀模式批量识别样品，支持断点续跑
- 比对与甲基化提取一体化，可控制是否忽略重叠 reads
- 自动生成样品级与项目级汇总报告及运行日志

## 快速开始 | Quick Start

```bash
biopytools bismark -g genome.fa -i ./raw_data -o ./bismark_results
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-g, --genome` | 参考基因组 FASTA 文件路径 |
| `-i, --input` | 原始 FASTQ 数据目录（程序在此目录内匹配 R1/R2） |
| `-o, --output-dir` | 主输出目录（不存在时自动创建） |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-p, --pattern` | `_1_clean.fq.gz` | R1 文件的后缀模式（R2 自动由 R1 推断） |
| `-t, --threads` | `12` | 比对/提取使用的线程数 |
| `--sort-buffer` | `400G` | 甲基化提取步骤的排序缓冲区大小（如 `400G`、`128G`） |
| `--include-overlap` | False | 包含 R1/R2 重叠读段（默认会通过 `--no_overlap` 排除） |
| `--version` | - | 显示版本信息 |

（运行 `biopytools bismark -h` 查看完整参数列表）

## 输出 | Output

输出目录典型结构：

- `bismark_index/`：Bismark 基因组索引（如不存在会先构建）
- `01_mapping/`：每样品的 `.bam` 比对结果
- `02_results/`：每样品的甲基化提取结果（含 CpG/CHG/CHH context 报告与 bedGraph）
- `03_summary/`：项目级汇总报告（比对率、甲基化率等）
- `99_logs/`：运行日志

## 依赖 | Dependencies

- `bismark`、`bismark_genome_preparation`、`bismark_methylation_extractor`（Bismark 套件）
- `bowtie2`（Bismark 比对后端）
- `samtools`（BAM 处理，Bismark 通常会随附）
- 推荐通过 `conda install -c bioconda bismark bowtie2 samtools` 安装

## 引用 | Citation

- Krueger F, Andrews SR. Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. *Bioinformatics*. 2011, 27(11):1571-1572.
- Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. *Nature Methods*. 2012, 9(4):357-359.

## 相关链接 | References

- [Bismark 官方文档](https://github.com/FelixKrueger/Bismark)
- [Bowtie2 主页](http://bowtie-bio.sourceforge.net/bowtie2/)
- [项目主页](https://github.com/lixiang117423/biopytools)
