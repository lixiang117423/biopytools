# BAM 文件批量统计分析 | BAM File Batch Statistics

**批量处理 BAM 文件，输出全局统计长表、染色体级别统计与基因组统计 JSON | Batch-process BAM files and produce summary TSV, per-chromosome TSV, and genome stats JSON.**

## 功能概述 | Overview

- 批量扫描单个 BAM 文件或目录下的所有 `.bam` 文件
- 多进程并行分析，每个 BAM 内部按多线程调用 samtools/bedtools
- 六大统计模块可独立开关：比对、覆盖度、序列特征、插入片段、重复、变异
- 按染色体拆分统计，并保存基因组级别指标 JSON
- MAPQ、插入片段长度、目标 BED 区域可自定义

## 快速开始 | Quick Start

```bash
# 分析目录下所有BAM
biopytools bam-stats -i ./bam_files -o result.summary.tsv

# 分析单个BAM并启用GC bias(需参考基因组)
biopytools bam-stats -i sample.bam -o sample.summary.tsv -g reference.fa
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | BAM 文件或包含 BAM 文件的目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `bam_stats.summary.tsv` | 全局统计输出文件 |
| `-t, --threads` | `12` | samtools 线程数 |
| `-p, --processes` | `16` | 并行处理的样本数 |
| `-g, --reference` | 无 | 参考基因组 FASTA (用于 GC bias 等) |
| `--bed-file` | 无 | 目标区域 BED 文件 |
| `--min-mapq` | `20` | 最小 MAPQ 阈值 |
| `--max-insert` | `1000` | 最大插入片段长度 |
| `--skip-alignment` | off | 跳过比对统计 |
| `--skip-coverage` | off | 跳过覆盖度统计 |
| `--skip-sequence` | off | 跳过序列特征统计 |
| `--skip-insert` | off | 跳过插入片段统计 |
| `--skip-duplicate` | off | 跳过重复统计 |
| `--skip-variation` | off | 跳过变异统计 |

（运行 `biopytools bam-stats -h` 查看完整参数列表）

## 输出 | Output

输出目录下生成以下文件（`{prefix}` 为 `--output` 文件主名）：

- `{prefix}.summary.tsv`：全局统计长表，每行一个样本，包含 Total_Reads、Map_Rate(%)、Avg_Depth、Cov_Bases_Rate(%) 等
- `{prefix}.per_chromosome.tsv`：染色体级别统计
- `{prefix}.genome_stats.json`：基因组级统计 JSON
- 运行日志保存于输出目录

## 依赖 | Dependencies

- samtools
- bedtools
- Python 包：pandas、tqdm

## 引用 | Citation

- Li H. et al. The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 2009.
- Quinlan A.R. & Hall I.M. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, 2010.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
