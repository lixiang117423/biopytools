# 🧬 K-mer提取工具 | K-mer Extraction Toolkit

## 📋 项目介绍 | Project Introduction

| 中文 | English |
|------|---------|
| K-mer提取工具是一个高性能的生物信息学工具，用于从FASTA和FASTQ文件中提取k-mer序列。该工具基于两个优秀的底层工具：对FASTA文件使用unikmer，对FASTQ文件使用jellyfish，提供了统一的接口和强大的功能。 | K-mer Extraction Toolkit is a high-performance bioinformatics tool for extracting k-mer sequences from FASTA and FASTQ files. This tool is based on two excellent underlying tools: unikmer for FASTA files and jellyfish for FASTQ files, providing a unified interface and powerful functionality. |

## ✨ 主要特性 | Key Features

| 功能 Feature | 中文描述 | English Description |
|--------------|----------|---------------------|
| 🔧 智能文件类型检测 | 自动识别FASTA/FASTQ文件格式，无需手动指定 | Automatic FASTA/FASTQ file format detection without manual specification |
| 🚀 高性能并行处理 | 支持多线程处理，充分利用计算资源 | Multi-threaded processing for optimal resource utilization |
| 🧬 双工具引擎 | FASTA文件使用unikmer，FASTQ文件使用jellyfish | Dual-engine approach: unikmer for FASTA, jellyfish for FASTQ |
| 🔗 配对文件支持 | 智能匹配双端测序文件（R1/R2） | Intelligent paired-end sequencing file matching (R1/R2) |
| 📁 批量处理 | 支持多文件同时处理和目录扫描 | Batch processing and directory scanning support |
| 🎯 灵活输出格式 | 支持FASTA、BED格式输出 | Flexible output formats: FASTA and BED |
| 🗜️ 压缩文件支持 | 原生支持.gz压缩文件 | Native support for .gz compressed files |
| 📊 详细统计信息 | 提供完整的k-mer提取统计报告 | Comprehensive k-mer extraction statistics |

## 🛠️ 安装要求 | Installation Requirements

### 依赖软件 | Required Software

| 软件 Software | 用途 Purpose | 安装命令 Installation |
|---------------|--------------|----------------------|
| **unikmer** | FASTA文件k-mer提取 \| FASTA k-mer extraction | `conda install -c bioconda unikmer` |
| **jellyfish** | FASTQ文件k-mer提取 \| FASTQ k-mer extraction | `conda install -c bioconda jellyfish` |
| **Python 3.7+** | 运行环境 \| Runtime environment | `python --version` |

### Python依赖 | Python Dependencies

```bash
# 无额外Python包依赖，仅使用标准库
# No additional Python packages required, uses standard library only
```

## 📦 安装步骤 | Installation Steps

### 方法1：直接使用 | Method 1: Direct Usage

```bash
# 1. 克隆或下载代码 | Clone or download code
git clone <repository-url>
cd kmer-extractor

# 2. 使用提取脚本提取模块文件 | Extract module files using extraction script
python extract_modules.py

# 3. 安装依赖软件 | Install required software
conda install -c bioconda unikmer jellyfish

# 4. 运行工具 | Run the tool
python run_kmer_extractor.py -i your_file.fasta -o results
```

### 方法2：包安装 | Method 2: Package Installation

```bash
# 如果打包为Python包 | If packaged as Python package
pip install kmer-extractor
```

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# FASTA文件处理 | FASTA file processing
python run_kmer_extractor.py -i genome.fasta -o results

# FASTQ文件处理 | FASTQ file processing  
python run_kmer_extractor.py -i reads.fastq.gz -o results -k 31
```

### Python API使用 | Python API Usage

```python
from kmer_extractor import KmerExtractor

# 创建提取器 | Create extractor
extractor = KmerExtractor(
    input_files=["sample.fasta"],
    output_dir="kmer_results",
    kmer_length=51,
    threads=16
)

# 运行提取 | Run extraction
extractor.run_extraction()
```

## 📖 详细用法 | Detailed Usage

### 命令行参数 | Command Line Arguments

| 参数 Parameter | 中文说明 | English Description | 默认值 Default |
|----------------|----------|---------------------|----------------|
| `-i, --input-files` | 输入文件路径（必需） | Input file paths (required) | - |
| `-o, --output-dir` | 输出目录 | Output directory | `./kmer_output` |
| `-k, --kmer-length` | K-mer长度(1-64) | K-mer length (1-64) | `51` |
| `-t, --threads` | 线程数 | Number of threads | `88` |
| `-m, --memory` | 内存限制(GB) | Memory limit (GB) | `880` |
| `--file-type` | 文件类型 (fasta/fastq) | File type (fasta/fastq) | 自动检测 auto-detect |
| `--fastq-pattern` | FASTQ匹配模式 | FASTQ matching pattern | - |
| `--output-bed` | 输出BED文件 | Output BED file | `False` |
| `--no-canonical` | 不使用canonical k-mer | Do not use canonical k-mers | `False` |
| `--no-compress` | 不压缩输出 | Do not compress output | `False` |
| `--no-keep-binary` | 不保留二进制文件 | Do not keep binary files | `False` |

### 使用示例 | Usage Examples

#### 1. 基本FASTA处理 | Basic FASTA Processing

```bash
# 中文：处理单个FASTA文件，使用51-mer
# English: Process single FASTA file with 51-mer
python run_kmer_extractor.py -i genome.fasta -o results
```

#### 2. FASTQ双端测序 | FASTQ Paired-end Sequencing

```bash
# 中文：处理双端测序FASTQ文件
# English: Process paired-end FASTQ files
python run_kmer_extractor.py \
    -i sample_R1.fastq.gz sample_R2.fastq.gz \
    -o results \
    --fastq-pattern "*_R1.fastq.gz" \
    -k 31 -t 16
```

#### 3. 批量文件处理 | Batch File Processing

```bash
# 中文：批量处理目录中的所有FASTA文件
# English: Batch process all FASTA files in directory
python run_kmer_extractor.py \
    -i /path/to/fasta/files/ \
    -o batch_results \
    -k 25 -t 32 \
    --output-bed
```

#### 4. 高性能配置 | High-performance Configuration

```bash
# 中文：高性能处理大型基因组
# English: High-performance processing for large genomes
python run_kmer_extractor.py \
    -i large_genome.fasta \
    -o hp_results \
    -k 51 -t 88 -m 880 \
    --output-bed
```

#### 5. 多样品FASTQ处理 | Multi-sample FASTQ Processing

```bash
# 中文：使用模式匹配处理多个样品
# English: Process multiple samples using pattern matching
python run_kmer_extractor.py \
    -i /data/samples/*.fq.gz \
    -o multi_sample_results \
    --fastq-pattern "*_1.clean.fq.gz" \
    -k 21 -t 24
```

## 📁 输出文件 | Output Files

### 输出文件类型 | Output File Types

| 文件类型 File Type | 说明 Description | 适用输入 Applicable Input |
|-------------------|------------------|---------------------------|
| `{base_name}.fasta` | FASTA格式的k-mer序列 \| K-mer sequences in FASTA format | 所有输入 All inputs |
| `{base_name}.bed` | BED格式的k-mer位置信息 \| K-mer position information in BED format | 仅FASTA输入 FASTA input only |
| `{base_name}.jf` | Jellyfish二进制文件 \| Jellyfish binary file | FASTQ输入 FASTQ input |
| `{base_name}.unik` | Unikmer二进制文件 \| Unikmer binary file | FASTA输入 FASTA input |
| `kmer_extraction.log` | 详细日志文件 \| Detailed log file | 所有输入 All inputs |

### 输出格式示例 | Output Format Examples

#### FASTA输出 | FASTA Output

```
>kmer_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG

>OV12_1_51  
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
```

#### BED输出 | BED Output

```
OV12    0       51      ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
OV12    1       52      TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
```

## ⚙️ 配置选项 | Configuration Options

### FASTQ模式匹配 | FASTQ Pattern Matching

| 模式示例 Pattern Example | 匹配文件 Matched Files | 说明 Description |
|-------------------------|------------------------|------------------|
| `*_1.fq.gz` | sample_1.fq.gz, sample_2.fq.gz | 双端测序标准命名 Standard paired-end naming |
| `*_R1.fastq` | data_R1.fastq, data_R2.fastq | Illumina标准命名 Illumina standard naming |
| `*.clean.fq` | sample.clean.fq | 单端测序 Single-end sequencing |

### 性能调优 | Performance Tuning

| 参数 Parameter | 推荐值 Recommended | 说明 Description |
|----------------|-------------------|------------------|
| 线程数 Threads | CPU核心数 CPU cores | 根据系统配置调整 Adjust based on system |
| 内存 Memory | 可用内存的80% 80% of available | 避免系统卡顿 Avoid system hang |
| K-mer长度 Length | 21-51 | 根据研究需求选择 Choose based on research needs |

## 🔧 故障排除 | Troubleshooting

### 常见问题 | Common Issues

| 问题 Issue | 中文解决方案 | English Solution |
|------------|-------------|------------------|
| 依赖软件未找到 | 确保unikmer/jellyfish在PATH中 | Ensure unikmer/jellyfish are in PATH |
| 内存不足 | 减少线程数或降低内存限制 | Reduce threads or memory limit |
| 文件权限错误 | 检查输入文件和输出目录权限 | Check input file and output directory permissions |
| K-mer长度错误 | 使用1-64范围内的值 | Use values within 1-64 range |

### 调试模式 | Debug Mode

```bash
# 查看详细日志 | View detailed logs
tail -f kmer_output/kmer_extraction.log

# 检查依赖 | Check dependencies
unikmer version
jellyfish --version
```

## 📊 性能基准 | Performance Benchmarks

### 测试数据 | Test Data

| 数据类型 Data Type | 文件大小 File Size | 处理时间 Processing Time | 内存使用 Memory Usage |
|-------------------|-------------------|-------------------------|----------------------|
| 小型基因组 Small genome | 10MB FASTA | ~30秒 30s | ~1GB |
| 中型基因组 Medium genome | 100MB FASTA | ~5分钟 5min | ~8GB |
| 大型基因组 Large genome | 1GB FASTA | ~30分钟 30min | ~50GB |
| FASTQ测序数据 FASTQ data | 500MB FASTQ | ~10分钟 10min | ~20GB |

## 🤝 贡献指南 | Contributing

### 开发环境 | Development Environment

```bash
# 克隆仓库 | Clone repository
git clone <repository-url>
cd kmer-extractor

# 创建开发分支 | Create development branch
git checkout -b feature/new-feature

# 运行测试 | Run tests
python -m pytest tests/
```

### 提交规范 | Commit Guidelines

| 类型 Type | 中文说明 | English Description |
|-----------|----------|---------------------|
| feat | 新功能 | New feature |
| fix | 错误修复 | Bug fix |
| docs | 文档更新 | Documentation update |
| style | 代码格式 | Code formatting |
| refactor | 重构 | Code refactoring |
| test | 测试相关 | Testing |

## 📄 许可证 | License

```
MIT License - 详见LICENSE文件
MIT License - See LICENSE file for details
```
