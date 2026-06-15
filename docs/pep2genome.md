# 🧬 Protein2Genome (Pep2Genome)

> 🎯 蛋白质到基因组比对分析工具 | Protein to Genome Alignment Analysis Tool

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/biopytools/biopytools)
[![Python](https://img.shields.io/badge/python-3.7+-green.svg)](https://www.python.org/downloads/)

## 📋 概述 | Overview

**Pep2Genome** 是一个基于 [Miniprot](https://github.com/lh3/miniprot) 的蛋白质到基因组比对分析工具，用于将蛋白质序列比对到参考基因组，生成PAF格式的比对结果，并支持导出GFF3和BED格式的基因注释文件。

**Pep2Genome** is a protein-to-genome alignment analysis tool based on [Miniprot](https://github.com/lh3/miniprot), which aligns protein sequences to a reference genome, generates PAF format alignment results, and supports exporting gene annotations in GFF3 and BED formats.

### ✨ 主要功能 | Key Features

- 🔍 **蛋白质比对** | **Protein Alignment**: 使用Miniprot进行高效的蛋白质到基因组比对
- 📊 **统计分析** | **Statistics Analysis**: 生成详细的比对统计报告，包括一致性、覆盖度等指标
- 📝 **格式导出** | **Format Export**: 支持导出GFF3和BED格式的基因注释文件
- ⚡ **多线程支持** | **Multi-threading**: 支持多线程加速比对过程
- 🎨 **质量评估** | **Quality Assessment**: 基于mapping quality和序列一致性评估比对质量

## 🚀 快速开始 | Quick Start

### 安装依赖 | Install Dependencies

```bash
# 使用conda安装Miniprot | Install Miniprot using conda
conda create -n miniprot_v.1.0.0 -c bioconda miniprot
conda activate miniprot_v.1.0.0
```

### 基本用法 | Basic Usage

```bash
# 基本比对分析 | Basic alignment analysis
biopytools pep2genome \
  --genome genome.fa \
  --protein protein.fa \
  -o results

# 指定线程数 | Specify thread count
biopytools pep2genome \
  --genome genome.fa \
  --protein protein.fa \
  -o results \
  -t 24

# 只导出特定格式 | Export specific formats only
biopytools pep2genome \
  --genome genome.fa \
  --protein protein.fa \
  -o results \
  --no-bed
```

## 📖 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | Parameter | 说明 | Description |
|------|-----------|------|-------------|
| `--genome` | Genome | 基因组FASTA文件 | Genome FASTA file |
| `--protein` | Protein | 蛋白质FASTA文件 | Protein FASTA file |
| `-o, --output-dir` | Output | 输出目录 | Output directory |

### 可选参数 | Optional Parameters

| 参数 | Parameter | 默认值 | Default | 说明 | Description |
|------|-----------|--------|---------|------|-------------|
| `-t, --threads` | Threads | `12` | `12` | 线程数 | Number of threads |

### 输出选项 | Output Options

| 参数 | Parameter | 说明 | Description |
|------|-----------|------|-------------|
| `--no-gff3` | No GFF3 | 不导出GFF3格式文件 | Do not export GFF3 format |
| `--no-bed` | No BED | 不导出BED格式文件 | Do not export BED format |
| `--no-statistics` | No Statistics | 不生成统计报告 | Do not generate statistics report |

## 📂 输出文件 | Output Files

分析完成后，输出目录包含以下文件：

After analysis completion, the output directory contains:

| 文件 | File | 说明 | Description |
|------|-------|------|-------------|
| `alignment.paf` | PAF | Miniprot比对结果 | Miniprot alignment results |
| `alignment_statistics.txt` | Statistics | 统计报告 | Statistics report |
| `alignment.gff3` | GFF3 | GFF3格式注释 | GFF3 format annotations |
| `alignment.bed` | BED | BED格式注释 | BED format annotations |

### PAF格式说明 | PAF Format Description

PAF (Pairwise mApping Format) 格式包含12列必选字段和多个可选tags：

PAF format includes 12 required columns and multiple optional tags:

```
列1:  query_name (查询序列名称)
列2:  query_length (查询序列长度)
列3:  query_start (查询起始位置)
列4:  query_end (查询终止位置)
列5:  strand (链方向 +/-)
列6:  target_name (目标序列名称)
列7:  target_length (目标序列长度)
列8:  target_start (目标起始位置)
列9:  target_end (目标终止位置)
列10: residue_matches (匹配残基数)
列11: alignment_block_length (比对块长度)
列12: mapping_quality (映射质量)
```

### 统计报告说明 | Statistics Report Description

统计报告包含以下内容：

The statistics report includes:

- **总体统计** | Overall Statistics: 总比对记录数、唯一查询/目标序列数、链方向分布
- **查询序列统计** | Query Statistics: 每个查询序列的比对数、多重比对统计
- **目标序列统计** | Target Statistics: 每个目标序列的比对数分布
- **比对质量统计** | Alignment Quality: Mapping quality、序列一致性分布
- **覆盖度统计** | Coverage Statistics: 查询序列和目标序列的覆盖度分布

## 🔧 环境变量 | Environment Variables

| 变量 | Variable | 说明 | Description |
|------|----------|------|-------------|
| `MINIPROT_PATH` | Miniprot Path | Miniprot工具路径 | Miniprot tool path |

## 📝 示例 | Examples

### 示例1: 基本比对 | Example 1: Basic Alignment

```bash
biopytools pep2genome \
  --genome /data/genome.fa \
  --protein /data/proteins.fa \
  -o ./pep2genome_results
```

### 示例2: 高性能比对 | Example 2: High-performance Alignment

```bash
biopytools pep2genome \
  --genome /data/genome.fa \
  --protein /data/proteins.fa \
  -o ./pep2genome_results \
  -t 48
```

### 示例3: 仅生成统计和GFF3 | Example 3: Statistics and GFF3 Only

```bash
biopytools pep2genome \
  --genome /data/genome.fa \
  --protein /data/proteins.fa \
  -o ./pep2genome_results \
  --no-bed
```

## ⚙️ 系统要求 | Requirements

- **Python**: >= 3.7
- **Miniprot**: >= 0.1 (可通过Bioconda安装 | Available via Bioconda)
- **内存** | Memory: 建议8GB以上 | Recommended 8GB+
- **磁盘空间** | Disk Space: 根据基因组大小而定 | Depends on genome size

## 🐛 故障排除 | Troubleshooting

### 问题1: Miniprot未找到 | Issue 1: Miniprot Not Found

```bash
# 解决方案：安装Miniprot | Solution: Install Miniprot
conda create -n miniprot_v.1.0.0 -c bioconda miniprot
conda activate miniprot_v.1.0.0
```

### 问题2: 比对时间过长 | Issue 2: Alignment Too Slow

```bash
# 解决方案：增加线程数 | Solution: Increase threads
biopytools pep2genome --genome genome.fa --protein protein.fa -o results -t 48
```

### 问题3: 内存不足 | Issue 3: Out of Memory

```bash
# 解决方案：减小蛋白质数据集或增加系统内存 | Solution: Reduce protein dataset or increase system memory
```

## 📚 参考文献 | References

- Li, H. (2023). Miniprot: protein-to-genome alignment in compact C. *bioRxiv*. [https://doi.org/10.1101/2023.04.26.538400](https://doi.org/10.1103/2023.04.26.538400)
- [Miniprot GitHub Repository](https://github.com/lh3/miniprot)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## 🤝 贡献 | Contributing

欢迎提交Issue和Pull Request！

Issues and Pull Requests are welcome!

---

**开发团队** | **Development Team**: Biopytools Development Team
**版本** | **Version**: 1.0.0
**最后更新** | **Last Updated**: 2025
