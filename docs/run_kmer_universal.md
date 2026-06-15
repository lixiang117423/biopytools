# 🧬 K-mer Universal Analysis Toolkit | 通用K-mer分析工具包

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/) [![License](https://img.shields.io/badge/license-MIT-green.svg)](https://claude.ai/chat/LICENSE) [![KMC](https://img.shields.io/badge/KMC-3.2.4+-orange.svg)](https://github.com/refresh-bio/KMC) [![Documentation](https://img.shields.io/badge/docs-available-brightgreen.svg)](https://claude.ai/chat/docs/)

------

## 🎯 简介 | Introduction

通用K-mer分析工具包是一个**高性能、灵活**的生物信息学工具，专为大规模FASTA和FASTQ文件的k-mer分析而设计。该工具支持智能文件格式识别、自动角色分配，并提供多种分析功能，特别适用于基因组学、转录组学和宏基因组学研究。

The K-mer Universal Analysis Toolkit is a **high-performance, flexible** bioinformatics tool designed for large-scale k-mer analysis of FASTA and FASTQ files. It supports intelligent file format detection, automatic role assignment, and provides various analysis capabilities, particularly suitable for genomics, transcriptomics, and metagenomics research.

### 🎨 设计理念 | Design Philosophy

- **🔄 双重角色**: 文件既可作为k-mer库来源，也可作为查询目标
- **🤖 智能识别**: 自动识别文件格式和最优角色分配
- **⚡ 高性能**: 基于KMC3的超高速k-mer处理
- **🎛️ 灵活配置**: 支持多种分析策略和输出格式

------

## ✨ 主要特性 | Key Features

### 🔍 智能文件处理 | Intelligent File Processing

- **格式识别**: 自动识别FASTA/FASTQ格式，支持.gz/.bz2/.xz压缩文件
- **角色分配**: 智能分配文件作为k-mer源或查询目标，支持5种分配策略
- **批量处理**: 支持通配符、目录扫描和文件列表

### ⚡ 高性能计算 | High-Performance Computing

- **KMC3集成**: 利用目前最快的k-mer计数工具
- **并行处理**: 支持88核心并行，1TB内存优化
- **大数据支持**: 可处理PB级数据（测试：200GB FASTA + 20TB FASTQ）
- **内存优化**: 智能内存管理，支持流式处理

### 📊 灵活分析 | Flexible Analysis

- **位置追踪**: FASTA文件保留完整位置信息（序列名、起始、终止位置）
- **样本管理**: FASTQ文件智能样本命名和计数
- **滑窗分析**: 可配置窗口大小的区域分析（默认500kb）
- **多种统计**: 丰度矩阵、存在/缺失矩阵、比例统计

### 📁 多样输出 | Multiple Output Formats

- **📄 FASTA格式**: K-mer库文件，标准序列格式
- **📊 CSV格式**: 丰度和存在矩阵，便于数据分析
- **📝 TXT格式**: 摘要报告和统计信息
- **🔄 JSON格式**: 结构化结果导出

------

## ⚙️ 安装 | Installation

### 📋 系统要求 | System Requirements

| 组件         | 最小要求   | 推荐配置                  |
| ------------ | ---------- | ------------------------- |
| **操作系统** | Linux/Unix | Ubuntu 20.04+ / CentOS 8+ |
| **Python**   | 3.8+       | 3.10+                     |
| **内存**     | 16GB       | 512GB-1TB                 |
| **CPU**      | 8核心      | 64-88核心                 |
| **存储**     | SSD推荐    | NVMe SSD                  |

### 🔧 依赖安装 | Dependencies Installation

#### 1️⃣ 安装KMC3 (必需) | Install KMC3 (Required)

```bash
# 📦 方法1: 使用conda (推荐)
conda install -c bioconda kmc

# 🔨 方法2: 从源码编译
git clone https://github.com/refresh-bio/KMC.git
cd KMC && make
export PATH=$PATH:$(pwd)/bin

# ✅ 验证安装
kmc -h
```

#### 2️⃣ 安装Python依赖 | Install Python Dependencies

```bash
# 📥 基础依赖
pip install biopython>=1.80 pandas>=1.5.0 numpy>=1.20.0 pyyaml>=6.0

# 📋 或使用requirements文件
pip install -r requirements.txt
```

### 📦 安装工具包 | Install Toolkit

```bash
# 🚀 方法1: 开发安装 (推荐)
git clone https://github.com/biopytools/kmer-universal.git
cd kmer-universal
pip install -e .

# 📥 方法2: PyPI安装
pip install biopytools-kmer-universal

# ✅ 验证安装
kmer-universal --version
```

------

## 🚀 快速开始 | Quick Start

### 🏃‍♂️ 最简单用法 | Simplest Usage

```bash
# 🔄 自动识别模式 - 让工具自动决定谁是k-mer源，谁是查询目标
kmer-universal --input genome.fa sample.fq --output results/
```

### 🎯 明确角色模式 | Explicit Role Mode

```bash
# 📋 明确指定：基因组作为k-mer源，样本作为查询目标
kmer-universal --kmer-sources genome.fa \
               --query-targets sample1.fq sample2.fq \
               --output results/
```

### ⚙️ 配置文件模式 | Configuration File Mode

创建配置文件 `analysis.yaml`:

```yaml
# 🧬 K-mer分析配置
kmer_size: 51
threads: 88
memory_gb: 880

# 📁 输入文件
input_paths:
  - "/data/genomes/*.fa"
  - "/data/samples/*.fq.gz"

# 🎛️ 分析参数
assignment_strategy: "intelligent"
window_sizes: [500000]
output_formats: ["fasta", "csv", "txt"]
output_dir: "results/"
# 🔧 使用配置文件运行
kmer-universal --config analysis.yaml
```

------

## 📝 使用示例 | Usage Examples

### 🧪 基础分析 | Basic Analysis

```bash
# 🔬 单个基因组vs单个样本
kmer-universal -i reference.fa sample.fastq -o basic_analysis/

# 🧬 指定k-mer大小和线程数
kmer-universal -i *.fa *.fq -o results/ -k 31 -t 64

# 📊 详细输出模式
kmer-universal -i *.fa *.fq -o results/ --verbose
```

### 🏭 批量处理 | Batch Processing

```bash
# 📂 处理整个目录
kmer-universal --input /data/genomes/ /data/samples/ \
               --assignment-strategy type_based \
               --output batch_results/

# 🌟 通配符批量处理
kmer-universal --input /project/*/genomes/*.fa /project/*/samples/*.fq.gz \
               --output large_scale_analysis/ \
               --threads 88 --memory 880
```

### 🎛️ 高级配置 | Advanced Configuration

```bash
# 🪟 自定义滑窗分析
kmer-universal -i ref.fa samples/*.fq \
               --window-sizes 100000 500000 1000000 \
               --kmer-size 51 \
               --output sliding_window_analysis/

# 💾 内存和缓存优化
kmer-universal -i large_dataset/ \
               --memory 1000 \
               --cache-size 200 \
               --chunk-size 5.0 \
               --threads 88 \
               --output optimized_results/

# 🤝 交互式角色分配
kmer-universal --input mixed_files/ \
               --interactive \
               --output interactive_results/
```

### 📊 不同角色分配策略 | Different Assignment Strategies

```bash
# 📏 基于文件大小分配（大文件作为查询目标）
kmer-universal --input mixed_data/ \
               --assignment-strategy size_based \
               --size-threshold 2.0 \
               --output size_based_results/

# 📁 基于文件类型分配（FASTA→源，FASTQ→目标）
kmer-universal --input mixed_data/ \
               --assignment-strategy type_based \
               --output type_based_results/

# 🧠 智能混合策略（推荐）
kmer-universal --input mixed_data/ \
               --assignment-strategy intelligent \
               --output intelligent_results/
```

------

## 📊 输出文件 | Output Files

分析完成后，在输出目录中生成以下文件：

After analysis completion, the following files are generated in the output directory:

### 📁 标准输出文件 | Standard Output Files

| 文件名                 | 格式  | 描述            | Description             |
| ---------------------- | ----- | --------------- | ----------------------- |
| `*_library.fasta`      | FASTA | 📄 K-mer库文件   | K-mer library file      |
| `*_abundance.csv`      | CSV   | 📊 丰度矩阵      | Abundance matrix        |
| `*_presence.csv`       | CSV   | ✅ 存在/缺失矩阵 | Presence/absence matrix |
| `*_sliding_window.csv` | CSV   | 🪟 滑窗分析结果  | Sliding window results  |
| `*_summary.txt`        | TXT   | 📋 分析摘要报告  | Analysis summary report |
| `kmer_analysis.log`    | LOG   | 📝 详细日志文件  | Detailed log file       |

### 📄 K-mer库文件格式 | K-mer Library Format

```fasta
>chr1_1000_1051               # FASTA来源: 序列名_起始位置_终止位置
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>sample1_kmer_000001          # FASTQ来源: 样本名_kmer_序号
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>chr2_5000_5051
TTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGCAATTGC
```

### 📊 丰度矩阵格式 | Abundance Matrix Format

**FASTA来源的k-mer**:

```csv
seq_name,start_pos,end_pos,kmer_seq,sample1,sample2,sample3
chr1,1000,1051,ATCGATCGATC...,5,0,12
chr1,1001,1052,TCGATCGATCG...,3,2,8
chr2,5000,5051,TTGCAATTGCA...,0,7,4
```

**FASTQ来源的k-mer**:

```csv
sample_name,kmer_id,kmer_seq,target1,target2,target3
ref_sample,kmer_000001,GCTAGCTAGCT...,3,0,5
ref_sample,kmer_000002,CTAGCTAGCTA...,1,4,2
```

### 🪟 滑窗分析格式 | Sliding Window Format

```csv
seq_name,window_start,window_end,window_size,sample_name,total_kmers,present_kmers,presence_ratio
chr1,1,500000,500000,sample1,1234,856,0.693
chr1,250001,750000,500000,sample1,1456,923,0.634
chr1,500001,1000000,500000,sample1,1123,789,0.703
```

### 📋 摘要报告内容 | Summary Report Content

```txt
=== K-mer Analysis Summary Report ===

📊 Basic Information:
  K-mer size: 51
  Total unique k-mers: 2,456,789
  Target samples: 3

🧬 K-mer Sources:
  From FASTA files: 1,234,567
  From FASTQ files: 1,222,222

📈 Sample Statistics:
  sample1:
    Present k-mers: 1,876,543/2,456,789 (76.4%)
    Total abundance: 15,678,901
  sample2:
    Present k-mers: 1,654,321/2,456,789 (67.3%)
    Total abundance: 12,345,678

⚙️ Configuration:
  Threads: 88
  Memory: 880GB
  Window sizes: [500000]
  Output directory: results/
```

------

## 🔧 参数说明 | Parameter Description

### 📥 输入输出参数 | Input/Output Parameters

| 参数              | 简写   | 默认值                  | 描述                       | Description                        |
| ----------------- | ------ | ----------------------- | -------------------------- | ---------------------------------- |
| `--input`         | `-i`   | -                       | 🔄 输入文件路径（自动模式） | Input file paths (auto mode)       |
| `--kmer-sources`  | `--ks` | -                       | 📄 K-mer源文件（明确模式）  | K-mer source files (explicit mode) |
| `--query-targets` | `--qt` | -                       | 🎯 查询目标文件（明确模式） | Query target files (explicit mode) |
| `--output`        | `-o`   | `kmer_analysis_results` | 📁 输出目录                 | Output directory                   |
| `--config`        | `-c`   | -                       | ⚙️ 配置文件（YAML格式）     | Configuration file (YAML format)   |

### 🧬 K-mer参数 | K-mer Parameters

| 参数             | 简写   | 默认值  | 范围       | 描述                |
| ---------------- | ------ | ------- | ---------- | ------------------- |
| `--kmer-size`    | `-k`   | `51`    | 1-256      | 🔢 K-mer长度         |
| `--min-count`    | `--ci` | `1`     | ≥0         | 📊 最小k-mer计数     |
| `--max-count`    | `--cx` | `1e9`   | >min-count | 📈 最大k-mer计数     |
| `--no-canonical` | `--nc` | `False` | -          | 🚫 禁用k-mer标准形式 |

### 🖥️ 系统资源参数 | System Resource Parameters

| 参数         | 简写    | 默认值                | 推荐值      | 描述           |
| ------------ | ------- | --------------------- | ----------- | -------------- |
| `--threads`  | `-t`    | `88`                  | CPU核心数   | 🔄 线程数       |
| `--memory`   | `-m`    | `880`                 | 80%可用内存 | 💾 内存限制(GB) |
| `--temp-dir` | `--tmp` | `/tmp/kmer_universal` | SSD路径     | 📂 临时目录     |

### 🧮 分析参数 | Analysis Parameters

| 参数                    | 简写   | 默认值        | 选项                                                   | 描述               |
| ----------------------- | ------ | ------------- | ------------------------------------------------------ | ------------------ |
| `--assignment-strategy` | `--as` | `intelligent` | explicit/size_based/type_based/intelligent/interactive | 🎛️ 角色分配策略     |
| `--window-sizes`        | `--ws` | `[500000]`    | >0                                                     | 🪟 滑窗大小(bp)     |
| `--size-threshold`      | `--st` | `1.0`         | >0                                                     | 📏 文件大小阈值(GB) |
| `--no-positions`        | `--np` | `False`       | -                                                      | 🚫 禁用位置追踪     |

### 📤 输出格式参数 | Output Format Parameters

| 参数                  | 简写    | 默认值          | 选项          | 描述              |
| --------------------- | ------- | --------------- | ------------- | ----------------- |
| `--output-formats`    | `--of`  | `fasta,csv,txt` | fasta/csv/txt | 📄 输出格式        |
| `--no-library`        | `--nl`  | `False`         | -             | 🚫 跳过k-mer库输出 |
| `--no-sliding-window` | `--nsw` | `False`         | -             | 🚫 跳过滑窗分析    |

------

## ⚡ 性能优化 | Performance Optimization

### 💻 硬件配置建议 | Hardware Configuration Recommendations

#### 🏠 小型分析 | Small Analysis (<10GB数据)

```bash
# 💻 配置建议
kmer-universal -i small_dataset/ \
               --threads 16 \
               --memory 128 \
               --chunk-size 1.0 \
               --output small_results/
```

- **内存**: 128-256GB
- **CPU**: 16-32核心
- **存储**: SATA SSD

#### 🏢 中型分析 | Medium Analysis (10-100GB数据)

```bash
# 🖥️ 配置建议
kmer-universal -i medium_dataset/ \
               --threads 64 \
               --memory 512 \
               --chunk-size 2.0 \
               --cache-size 100 \
               --output medium_results/
```

- **内存**: 512GB-1TB
- **CPU**: 32-64核心
- **存储**: NVMe SSD

#### 🏭 大型分析 | Large Analysis (>100GB数据)

```bash
# 🚀 高性能配置
kmer-universal -i large_dataset/ \
               --threads 88 \
               --memory 880 \
               --chunk-size 5.0 \
               --cache-size 200 \
               --ram-only \
               --output large_results/
```

- **内存**: 1TB+
- **CPU**: 64-88核心
- **存储**: 多NVMe SSD RAID

### ⚙️ 性能调优参数 | Performance Tuning Parameters

#### 💾 内存优化 | Memory Optimization

```bash
# 🧠 智能内存管理
kmer-universal -i dataset/ \
               --memory 1000 \           # 总内存限制
               --cache-size 300 \        # 缓存大小
               --chunk-size 5.0 \        # 分片大小
               --strict-memory           # 严格内存控制
```

#### 💨 I/O优化 | I/O Optimization

```bash
# 💽 高速I/O配置
kmer-universal -i dataset/ \
               --temp-dir /nvme/tmp/ \   # 使用NVMe存储
               --enable-compression \    # 启用压缩
               --signature-length 11     # 增加签名长度
```

#### 🔄 并行优化 | Parallel Optimization

```bash
# ⚡ 最大并行化
kmer-universal -i dataset/ \
               --threads 88 \            # 最大线程数
               --ram-only \              # RAM模式
               --no-strict-memory        # 允许内存超限
```

### 📊 性能基准 | Performance Benchmarks

| 数据规模                 | 硬件配置   | 预计时间 | 内存使用 | 存储需求 |
| ------------------------ | ---------- | -------- | -------- | -------- |
| 1GB FASTA + 10GB FASTQ   | 32核/256GB | 30分钟   | 150GB    | 50GB     |
| 10GB FASTA + 100GB FASTQ | 64核/512GB | 2小时    | 400GB    | 300GB    |
| 200GB FASTA + 20TB FASTQ | 88核/1TB   | 8-15小时 | 900GB    | 5TB      |

------

## 🔍 故障排除 | Troubleshooting

### ❗ 常见错误 | Common Errors

#### 🚫 KMC未找到 | KMC Not Found

```bash
Error: KMC not found in PATH. Please install KMC3
```

**解决方案 | Solution**:

```bash
# ✅ 检查KMC安装
which kmc
kmc -h

# 🔧 重新安装KMC
conda install -c bioconda kmc
# 或 or
export PATH=$PATH:/path/to/kmc/bin
```

#### 💾 内存不足 | Out of Memory

```bash
Error: Cannot allocate memory
```

**解决方案 | Solution**:

```bash
# 🔽 减少内存使用
kmer-universal -i dataset/ \
               --memory 500 \           # 减少内存限制
               --chunk-size 1.0 \       # 减小分片大小
               --no-ram-only            # 使用磁盘模式

# 🧹 清理系统内存
sudo sync && sudo sysctl vm.drop_caches=3
```

#### 📁 文件权限错误 | File Permission Error

```bash
Error: Permission denied
```

**解决方案 | Solution**:

```bash
# 🔐 检查文件权限
ls -la input_files/
chmod 644 input_files/*

# 📂 检查输出目录权限
mkdir -p results/
chmod 755 results/
```

### 🧪 调试模式 | Debug Mode

```bash
# 🔍 干运行模式 - 查看将要执行的操作
kmer-universal -i files/ --dry-run

# 📝 详细输出模式
kmer-universal -i files/ --verbose

# 🔧 保留中间文件用于调试
kmer-universal -i files/ --keep-intermediate

# 📊 检查文件识别结果
kmer-universal -i mixed_files/ --interactive
```

### 🚨 紧急处理 | Emergency Procedures

#### ⏹️ 优雅停止分析 | Graceful Stop

```bash
# Ctrl+C 或发送SIGTERM信号
kill -TERM <process_id>
```

#### 🔄 恢复中断的分析 | Resume Interrupted Analysis

```bash
# 🔄 断点续传功能
kmer-universal -i dataset/ --resume --output same_output_dir/
```

#### 🧹 清理临时文件 | Cleanup Temporary Files

```bash
# 🗑️ 手动清理
rm -rf /tmp/kmer_universal/*
rm -rf output_dir/*.kmc_*
```

### 📞 获取帮助 | Getting Help

```bash
# 📋 查看完整帮助
kmer-universal --help

# 📖 查看版本信息
kmer-universal --version

# 💾 导出当前配置
kmer-universal --save-config current_config.yaml --dry-run
```

------

## 👥 开发指南 | Development Guide

### 🏗️ 项目结构 | Project Structure

```
biopytools/kmer_universal/
├── 📁 __init__.py              # 包初始化 | Package initialization
├── 📁 config.py                # 配置管理 | Configuration management
├── 📁 file_manager.py          # 文件管理和识别 | File management and detection
├── 📁 kmc_interface.py         # KMC3接口封装 | KMC3 interface wrapper
├── 📁 position_tracker.py      # 位置信息追踪 | Position tracking
├── 📁 output_manager.py        # 输出管理 | Output management
├── 📁 analyzer.py              # 主分析逻辑 | Main analysis logic
└── 📁 main.py                  # 命令行接口 | Command line interface
```

### 🔌 API使用 | API Usage

```python
from biopytools.kmer_universal import KmerAnalyzer, KmerConfig

# 🛠️ 创建自定义配置
config = KmerConfig(
    kmer_size=31,
    threads=64,
    memory_gb=512,
    assignment_strategy=AssignmentStrategy.INTELLIGENT
)

# 🚀 创建分析器
analyzer = KmerAnalyzer(config)

# 🔄 执行自动分析
results = analyzer.auto_analyze(
    input_paths=["/data/genomes/", "/data/samples/"],
    output_dir="custom_results/"
)

# 📊 获取结果信息
print(f"K-mer library size: {results['kmer_library_size']}")
print(f"Target samples: {results['target_samples']}")
```

### 🧪 扩展开发 | Extension Development

#### 📝 自定义分析插件 | Custom Analysis Plugin

```python
from biopytools.kmer_universal.analyzer import KmerAnalyzer

class CustomKmerAnalyzer(KmerAnalyzer):
    def custom_analysis(self, kmer_library, target_abundances):
        """🔬 自定义分析逻辑"""
        # 实现你的分析逻辑
        pass
```

#### 🔧 自定义输出格式 | Custom Output Format

```python
from biopytools.kmer_universal.output_manager import OutputManager

class CustomOutputManager(OutputManager):
    def write_custom_format(self, data, output_file):
        """📄 自定义输出格式"""
        # 实现你的输出逻辑
        pass
```

### 🧪 测试框架 | Testing Framework

```bash
# 🔬 运行单元测试
python -m pytest tests/ -v

# 📊 运行覆盖率测试
python -m pytest tests/ --cov=biopytools.kmer_universal

# ⚡ 运行性能测试
python tests/performance_tests.py
```

### 📚 文档生成 | Documentation Generation

```bash
# 📖 生成API文档
sphinx-build -b html docs/ docs/_build/

# 🔄 自动文档更新
sphinx-autobuild docs/ docs/_build/
```

------

## 🤝 贡献 | Contributing

我们欢迎所有形式的贡献！🎉

We welcome all forms of contributions! 🎉

### 🛠️ 贡献类型 | Types of Contributions

- 🐛 **Bug报告** | Bug Reports: 发现和报告问题
- 💡 **功能建议** | Feature Requests: 提出新功能想法
- 🔧 **代码贡献** | Code Contributions: 提交代码改进
- 📖 **文档改进** | Documentation: 改进文档和示例
- 🧪 **测试用例** | Test Cases: 添加测试用例

### 📋 贡献流程 | Contribution Process

1. **🍴 Fork项目仓库** | Fork the repository

2. 🌿 创建特性分支

    | Create a feature branch:

   ```bash
   git checkout -b feature/amazing-feature
   ```

3. ✍️ 提交更改

    | Commit your changes:

   ```bash
   git commit -m "✨ Add amazing feature"
   ```

4. 🚀 推送分支

    | Push to the branch:

   ```bash
   git push origin feature/amazing-feature
   ```

5. **📝 提交Pull Request** | Open a Pull Request

### 📏 代码规范 | Code Standards

```bash
# 🖤 代码格式化
black biopytools/

# 🔍 代码检查
flake8 biopytools/

# 🔬 类型检查
mypy biopytools/
```

### 🧪 开发环境设置 | Development Environment Setup

```bash
# 🔧 克隆开发版本
git clone https://github.com/your-username/kmer-universal.git
cd kmer-universal

# 🐍 创建虚拟环境
python -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\Scripts\activate   # Windows

# 📦 安装开发依赖
pip install -e ".[dev]"

# ✅ 运行测试确保环境正确
python -m pytest
```

------

## 📞 支持与社区 | Support & Community

### 💬 获取帮助 | Getting Help

- 📧 **Email**: biopytools@example.com
- 🐛 **Issues**: [GitHub Issues](https://github.com/biopytools/kmer-universal/issues)
- 📚 **文档**: [Read the Docs](https://biopytools-kmer-universal.readthedocs.io/)
- 💬 **讨论**: [GitHub Discussions](https://github.com/biopytools/kmer-universal/discussions)

------

## 📖 引用 | Citation

如果在研究中使用了此工具，请引用：

If you use this tool in your research, please cite:

```bibtex
@software{kmer_universal_2024,
  title={K-mer Universal Analysis Toolkit: High-Performance Flexible K-mer Analysis for Genomics},
  author={BioPyTools Development Team},
  year={2024},
  url={https://github.com/biopytools/kmer-universal},
  version={1.0.0},
  doi={10.5281/zenodo.1234567}
}
```

------

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](https://claude.ai/chat/LICENSE) 文件。

This project is licensed under the MIT License - see the [LICENSE](https://claude.ai/chat/LICENSE) file for details.

```
MIT License

Copyright (c) 2024 BioPyTools Development Team

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

------

## 🎯 路线图 | Roadmap

### 🚀 v1.0.0 (当前版本 | Current)

- ✅ 基础k-mer分析功能
- ✅ 智能文件识别
- ✅ 多种角色分配策略
- ✅ 滑窗分析

### 🔄 v1.1.0 (计划中 | Planned)

- 🔄 GPU加速支持
- 🔄 分布式计算支持
- 🔄 更多输出格式
- 🔄 可视化界面

### 🌟 v2.0.0 (未来 | Future)

- 🌟 机器学习集成
- 🌟 云平台支持
- 🌟 实时分析功能
- 🌟 Web界面

------

<div align="center">

## 🎉 感谢使用 K-mer Universal Analysis Toolkit! | Thank you for using K-mer Universal Analysis Toolkit!

**⭐ 如果这个项目对您有帮助，请给我们一个星标！**

**⭐ If this project helps you, please give us a star!**

------

[![GitHub stars](https://img.shields.io/github/stars/biopytools/kmer-universal.svg?style=social&label=Star)](https://github.com/biopytools/kmer-universal) [![GitHub forks](https://img.shields.io/github/forks/biopytools/kmer-universal.svg?style=social&label=Fork)](https://github.com/biopytools/kmer-universal/fork)

</div>