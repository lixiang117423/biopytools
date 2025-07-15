# Biohelpers

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![PyPI Version](https://img.shields.io/badge/pypi-0.2.0-orange.svg)](#)

**生物信息学常用工具包** | **Bioinformatics Common Toolkits**

一个集成了多种生物信息学分析工具的Python包，提供统一的API和命令行接口。

A Python package integrating various bioinformatics analysis tools with unified API and command-line interfaces.

---

## 安装 | Installation

### 从源码安装 | Install from source

```bash
git clone https://github.com/lixiang117423/biohelpers.git
cd biohelpers
pip install -e .
```

### 安装开发依赖 | Install development dependencies

```bash
pip install -e ".[dev]"
```

---

# fastp 质控模块

## 功能概述

fastp 质控模块是 biopytools 工具包中的高效 FASTQ 数据质量控制工具，支持单端和双端测序数据的批量处理。该模块封装了 fastp 工具，提供了便捷的批处理功能和灵活的参数配置。

### 主要特性

- 🚀 **高效批处理**：自动识别和处理整个目录下的 FASTQ 文件
- 🔧 **灵活配置**：支持多种质控参数的自定义设置
- 📊 **质量报告**：自动生成 HTML 和 JSON 格式的质控报告
- 💾 **双端支持**：完整支持双端测序数据（Paired-end）处理
- 🎯 **智能识别**：自动识别文件配对关系

## 安装方法

### 系统依赖

确保系统已安装 fastp：

```bash
# Ubuntu/Debian
sudo apt-get install fastp

# CentOS/RHEL
sudo yum install fastp

# macOS (使用 Homebrew)
brew install fastp

# 或者使用 conda
conda install -c bioconda fastp
```

### 安装 biopytools

```bash
# 克隆项目
git clone https://github.com/yourusername/biopytools.git
cd biopytools

# 安装包
pip install -e .

# 验证安装
run_fastp --help
```

## 使用方法

### 基本语法

```bash
run_fastp -i INPUT_DIR -o OUTPUT_DIR [OPTIONS]
```

### 必需参数

| 参数 | 说明 |
|------|------|
| `-i, --input-dir` | 输入原始 FASTQ 数据目录 |
| `-o, --output-dir` | 输出清洁 FASTQ 数据目录 |

### 可选参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--fastp-path` | `fastp` | fastp 可执行文件路径 |
| `-t, --threads` | `12` | 线程数 |
| `-q, --quality-threshold` | `30` | 质量阈值 |
| `-l, --min-length` | `50` | 最小长度 |
| `-u, --unqualified-percent` | `40` | 不合格碱基百分比阈值 |
| `-n, --n-base-limit` | `10` | N 碱基数量限制 |
| `--read1-suffix` | `_1.fq.gz` | Read1 文件后缀 |
| `--read2-suffix` | `_2.fq.gz` | Read2 文件后缀 |

## 使用示例

### 1. 基本用法

```bash
# 处理目录下的所有 FASTQ 文件
run_fastp -i ./raw_data -o ./clean_data
```

### 2. 自定义参数

```bash
# 使用更严格的质控标准
run_fastp -i ./raw_data -o ./clean_data \
    -q 35 \
    -l 75 \
    -u 30 \
    -t 16
```

### 3. 不同文件后缀

```bash
# 处理以 .R1.fastq.gz 和 .R2.fastq.gz 结尾的文件
run_fastp -i ./raw_data -o ./clean_data \
    --read1-suffix .R1.fastq.gz \
    --read2-suffix .R2.fastq.gz
```

### 4. 指定 fastp 路径

```bash
# 使用自定义路径的 fastp
run_fastp -i ./raw_data -o ./clean_data \
    --fastp-path /usr/local/bin/fastp
```

### 5. 完整参数示例

```bash
run_fastp \
    -i /path/to/raw_data \
    -o /path/to/clean_data \
    --fastp-path /usr/local/bin/fastp \
    -t 20 \
    -q 25 \
    -l 40 \
    -u 50 \
    -n 5 \
    --read1-suffix _R1.fq.gz \
    --read2-suffix _R2.fq.gz
```

## 输入文件格式

### 目录结构要求

```
raw_data/
├── sample1_1.fq.gz    # Read1
├── sample1_2.fq.gz    # Read2  
├── sample2_1.fq.gz
├── sample2_2.fq.gz
└── ...
```

### 支持的文件格式

- **.fq.gz** / **.fastq.gz**：压缩的 FASTQ 文件
- **.fq** / **.fastq**：未压缩的 FASTQ 文件

### 文件命名规则

程序会根据 `--read1-suffix` 和 `--read2-suffix` 参数自动识别配对文件：

- `sample_1.fq.gz` ↔ `sample_2.fq.gz`
- `sample_R1.fastq.gz` ↔ `sample_R2.fastq.gz`
- `sample.R1.fq.gz` ↔ `sample.R2.fq.gz`

## 输出结果

### 输出目录结构

```
clean_data/
├── sample1_1.clean.fq.gz     # 清洁的 Read1 文件
├── sample1_2.clean.fq.gz     # 清洁的 Read2 文件
├── sample1.fastp.html        # HTML 质控报告
├── sample1.fastp.json        # JSON 质控报告
├── sample2_1.clean.fq.gz
├── sample2_2.clean.fq.gz
├── sample2.fastp.html
├── sample2.fastp.json
└── batch_summary.txt         # 批处理总结报告
```

### 质控报告说明

- **HTML 报告**：可视化的质控结果，包含质量分布图、GC含量等
- **JSON 报告**：机器可读的质控统计数据
- **批处理总结**：所有样本的处理状态和统计信息

## 质控参数说明

### 质量阈值 (-q, --quality-threshold)

- **默认值**：30
- **说明**：Phred 质量值阈值，低于此值的碱基被认为是低质量碱基
- **建议值**：
  - 严格：35+
  - 标准：30
  - 宽松：20-25

### 最小长度 (-l, --min-length)

- **默认值**：50
- **说明**：过滤后序列的最小长度，短于此长度的序列将被丢弃
- **建议值**：
  - RNA-seq：50-75
  - DNA-seq：30-50
  - 16S rRNA：200+

### 不合格碱基百分比 (-u, --unqualified-percent)

- **默认值**：40
- **说明**：如果序列中低质量碱基的百分比超过此阈值，整条序列将被丢弃
- **建议值**：30-50%

### N 碱基限制 (-n, --n-base-limit)

- **默认值**：10
- **说明**：序列中允许的最大 N 碱基数量
- **建议值**：5-15

## 性能优化

### 线程设置

```bash
# 根据 CPU 核心数设置线程
run_fastp -i input -o output -t $(nproc)

# 或者设置为核心数的 80%
run_fastp -i input -o output -t $(($(nproc) * 4 / 5))
```

### 内存使用

- 每个线程大约使用 500MB-1GB 内存
- 建议总内存使用量不超过系统内存的 80%

## 故障排除

### 常见问题

**1. 找不到 fastp 命令**
```bash
# 检查 fastp 是否安装
which fastp

# 如果未安装，使用 conda 安装
conda install -c bioconda fastp
```

**2. 权限错误**
```bash
# 确保输出目录有写入权限
chmod 755 /path/to/output_dir
```

**3. 文件未找到**
```bash
# 检查输入目录是否存在
ls -la /path/to/input_dir

# 检查文件后缀是否正确
ls /path/to/input_dir/*_1.fq.gz
```

**4. 内存不足**
```bash
# 减少线程数
run_fastp -i input -o output -t 4
```

### 调试模式

```bash
# 查看详细输出
run_fastp -i input -o output --verbose

# 检查 fastp 版本
fastp --version
```

## 最佳实践

### 1. 质控前检查

```bash
# 检查原始数据质量
fastqc raw_data/*.fq.gz -o qc_reports/

# 统计文件数量
find raw_data -name "*_1.fq.gz" | wc -l
```

### 2. 参数选择建议

| 应用场景 | 质量阈值 | 最小长度 | 不合格碱基% |
|----------|----------|----------|-------------|
| RNA-seq | 25-30 | 50-75 | 40-50 |
| WGS | 30-35 | 50-100 | 30-40 |
| 16S rRNA | 25 | 200+ | 50 |
| ChIP-seq | 20-25 | 30-50 | 50 |

### 3. 质控后验证

```bash
# 检查处理结果
fastqc clean_data/*.clean.fq.gz -o qc_reports_after/

# 比较处理前后的统计
multiqc qc_reports/ qc_reports_after/
```

## Python API 使用

```python
from biopytools.fastp import FastpProcessor, FastpConfig

# 创建配置
config = FastpConfig(
    input_dir="./raw_data",
    output_dir="./clean_data",
    threads=16,
    quality_threshold=30,
    min_length=50
)

# 运行处理
processor = FastpProcessor(config)
results = processor.run_batch_processing()

# 查看结果
print(f"处理了 {results['total_samples']} 个样本")
print(f"成功：{results['success_count']}")
print(f"失败：{results['failed_count']}")
```

## 更新日志

### v1.0.0
- 初始版本发布
- 支持双端测序数据批量处理
- 集成 fastp 质控功能
- 生成 HTML 和 JSON 报告

---

## 相关链接

- [fastp GitHub](https://github.com/OpenGene/fastp)
- [fastp 文档](https://github.com/OpenGene/fastp#usage)
- [FASTQ 格式说明](https://en.wikipedia.org/wiki/FASTQ_format)

---

## 贡献指南 | Contributing

1. Fork本仓库
2. 创建功能分支 (`git checkout -b feature/amazing-feature`)
3. 提交更改 (`git commit -m 'Add amazing feature'`)
4. 推送到分支 (`git push origin feature/amazing-feature`)
5. 创建Pull Request

---

## 许可证 | License

本项目使用MIT许可证 - 查看 [LICENSE](LICENSE) 文件了解详情

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 支持与反馈 | Support & Feedback

- 📧 Email: biohelpers@example.com
- 🐛 Issues: [GitHub Issues](https://github.com/yourusername/biohelpers/issues)
- 📖 Documentation: [Read the Docs](https://biohelpers.readthedocs.io)

---

# 更新日志 | Changelog

## v1.1.0 (2025-07-15)
- 增加run-fastp模块
