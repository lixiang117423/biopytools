# FASTP 质控分析模块

**专业的FASTQ数据质量控制工具 | Professional FASTQ Data Quality Control Tool**

## 功能概述 | Overview

FASTP 质控分析模块是一个高效的FASTQ数据质量控制工具，基于fastp软件构建，提供从原始数据质量控制到清洁数据输出的完整流程。支持单末端和双末端测序数据的自动化批量处理、灵活的质量控制参数配置、详细的质控报告生成，适用于各种高通量测序数据的质量控制分析。

## 主要特性 | Key Features

- **高效批处理**: 自动识别和处理整个目录下的FASTQ文件，支持单文件和目录模式
- **灵活配置**: 支持多种质控参数自定义设置(质量阈值、长度过滤、N碱基限制等)
- **智能配对**: 自动识别Read1和Read2配对文件，支持多种文件命名模式
- **详细报告**: 自动生成HTML和JSON格式的可视化质控报告
- **双端支持**: 完整支持单末端和双末端测序数据处理
- **质量控制**: 多维度质量过滤(质量值、序列长度、N碱基、不合格碱基百分比)
- **日志记录**: 完整的处理过程日志和错误追踪
- **高效处理**: 优化的多线程处理流程，支持大规模测序数据

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 处理整个目录的双末端测序数据
biopytools fastp \
    -i raw_data/ \
    -o clean_data/

# 使用自定义质控参数
biopytools fastp \
    -i raw_data/ \
    -o clean_data/ \
    -q 35 \
    -l 75 \
    -t 16
```

### 高级用法 | Advanced Usage

```bash
# 处理单末端测序数据
biopytools fastp \
    -i single_end_data/ \
    -o clean_data/ \
    --single-end

# 使用自定义文件后缀和质控参数
biopytools fastp \
    -i raw_data/ \
    -o clean_data/ \
    --read1-suffix .R1.fastq.gz \
    --read2-suffix .R2.fastq.gz \
    -q 25 \
    -l 50 \
    -u 30 \
    -n 5 \
    -t 24

# 强制覆盖已存在文件
biopytools fastp \
    -i raw_data/ \
    -o clean_data/ \
    --force
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入原始FASTQ数据目录或文件 | `-i raw_data/` 或 `-i sample_1.fq.gz` |
| `-o, --output-dir` | 输出清洁FASTQ数据目录 | `-o clean_data/` |

### 软件配置 | Software Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--fastp-path` | `fastp` | fastp可执行文件路径 |

### 质控参数 | Quality Control Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `12` | 线程数 |
| `-q, --quality-threshold` | `30` | 质量阈值(Phred质量值，0-50) |
| `-l, --min-length` | `50` | 最小长度要求 |
| `-u, --unqualified-percent` | `40` | 不合格碱基百分比阈值(0-100) |
| `-n, --n-base-limit` | `10` | N碱基数量限制 |

### 文件模式 | File Patterns

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--read1-suffix` | `_1.fq.gz` | Read1文件后缀(单末端模式也使用此参数) |
| `--read2-suffix` | `_2.fq.gz` | Read2文件后缀 |
| `--single-end` | `False` | 启用单末端模式 |

### 日志选项 | Logging Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-v, --verbose` | `0` | 详细输出模式(-v: INFO, -vv: DEBUG) |
| `--quiet` | `False` | 静默模式(仅输出ERROR) |
| `--log-level` | `INFO` | 日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL) |
| `--log-file` | `None` | 日志文件路径 |

### 执行选项 | Execution Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-f, --force` | `False` | 强制覆盖已存在文件 |
| `--dry-run` | `False` | 模拟运行(不实际执行) |

## 输入文件格式 | Input File Formats

### 目录结构要求 | Directory Structure Requirements

**双末端测序数据 | Paired-end Data:**

```
raw_data/
├── sample1_1.fq.gz    # Read1文件
├── sample1_2.fq.gz    # Read2配对文件
├── sample2_1.fq.gz
├── sample2_2.fq.gz
└── ...
```

**单末端测序数据 | Single-end Data:**

```
raw_data/
├── sample1.fq.gz
├── sample2.fq.gz
└── ...
```

### 支持的文件格式 | Supported File Formats

- **压缩格式**: `.fq.gz`, `.fastq.gz`, `.fq.gzip`
- **未压缩格式**: `.fq`, `.fastq`

### 文件命名规则 | File Naming Rules

程序自动识别多种配对文件命名模式：

| Read1模式 | Read2模式 | 样本名 |
|-----------|-----------|--------|
| `sample_1.fq.gz` | `sample_2.fq.gz` | `sample` |
| `sample.R1.fastq.gz` | `sample.R2.fastq.gz` | `sample` |
| `sample_R1.fq.gz` | `sample_R2.fq.gz` | `sample` |
| `sample.r1.fq.gz` | `sample.r2.fq.gz` | `sample` |
| `sample.read1.fq.gz` | `sample.read2.fq.gz` | `sample` |

## 使用示例 | Usage Examples

### 示例1：基本质控流程 | Example 1: Basic Quality Control Pipeline

```bash
# 处理标准的双末端测序数据
biopytools fastp \
    -i /data/raw sequencing/ \
    -o /data/clean_data/
```

### 示例2：严格质控标准 | Example 2: Strict Quality Control Standards

```bash
# 使用更严格的质控参数
biopytools fastp \
    -i raw_data/ \
    -o clean_data/ \
    -q 35 \
    -l 75 \
    -u 30 \
    -n 5 \
    -t 24
```

### 示例3：不同文件命名模式 | Example 3: Different File Naming Patterns

```bash
# 处理以.R1.fastq.gz和.R2.fastq.gz结尾的文件
biopytools fastp \
    -i raw_data/ \
    -o clean_data/ \
    --read1-suffix .R1.fastq.gz \
    --read2-suffix .R2.fastq.gz
```

### 示例4：单末端测序数据 | Example 4: Single-end Sequencing Data

```bash
# 处理单末端测序数据
biopytools fastp \
    -i single_end_raw/ \
    -o single_end_clean/ \
    --single-end \
    --read1-suffix .fq.gz
```

### 示例5：处理单个文件 | Example 5: Process Single File

```bash
# 处理单个FASTQ文件(自动检测配对文件)
biopytools fastp \
    -i sample_1.fq.gz \
    -o clean_data/
```

### 示例6：强制覆盖和调试 | Example 6: Force Overwrite and Debugging

```bash
# 强制覆盖已存在文件并启用详细日志
biopytools fastp \
    -i raw_data/ \
    -o clean_data/ \
    --force \
    -vv
```

### 示例7：指定fastp路径 | Example 7: Specify Custom fastp Path

```bash
# 使用自定义安装的fastp
biopytools fastp \
    -i raw_data/ \
    -o clean_data/ \
    --fastp-path /usr/local/bin/fastp \
    -t 20
```

### 示例8：模拟运行测试 | Example 8: Dry Run Test

```bash
# 测试命令但不实际执行
biopytools fastp \
    -i raw_data/ \
    -o clean_data/ \
    --dry-run \
    -v
```

## 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
clean_data/
├── sample1_1.clean.fq.gz           # 质控后的Read1文件
├── sample1_2.clean.fq.gz           # 质控后的Read2文件
├── sample2_1.clean.fq.gz
├── sample2_2.clean.fq.gz
├── fastp_processing_summary.txt    # 批处理总结报告
└── fastp_reports/                  # 质控报告目录
    ├── sample1.html                # HTML质控报告
    ├── sample1.json                # JSON质控数据
    ├── sample2.html
    └── sample2.json
```

### 质控报告说明 | Quality Control Report Description

#### HTML报告 | HTML Report
- **质量分布图**: 读取质量分布随位置变化
- **GC含量图**: GC含量分布统计
- **序列长度分布**: 过滤前后序列长度对比
- **接头含量**: 接头序列检测和统计
- **过读序列**: 低质量碱基位置统计

#### JSON报告 | JSON Report
- 机器可读的质控统计数据
- 便于后续自动化分析
- 包含所有质控指标的数值

#### 总结报告 | Summary Report
- 总样本数和成功率
- 处理参数记录
- 输入输出路径信息
- 处理时间统计

## 质控参数详解 | Quality Control Parameters Explained

### 质量阈值 (-q, --quality-threshold)

- **参数含义**: Phred质量分数阈值，低于此值的碱基被认为是低质量碱基
- **取值范围**: 0-50
- **默认值**: 30

| 应用场景 | 推荐值 | 说明 |
|----------|--------|------|
| 严格质控 | 35+ | 适用于高质量数据，要求更高的准确性 |
| 标准质控 | 30 | 适用于大多数测序数据 |
| 宽松质控 | 20-25 | 适用于质量较差或珍贵的数据 |

### 最小长度 (-l, --min-length)

- **参数含义**: 过滤后序列的最小长度，短于此长度的序列将被丢弃
- **默认值**: 50

| 测序类型 | 推荐长度 | 说明 |
|----------|----------|------|
| RNA-seq | 50-75 | 转录本测序通常较短 |
| WGS | 50-100 | 全基因组测序 |
| 16S rRNA | 200+ | 16S rRNA全长扩增 |
| ChIP-seq | 30-50 | 染色质免疫沉淀测序 |

### 不合格碱基百分比 (-u, --unqualified-percent)

- **参数含义**: 如果序列中低质量碱基的百分比超过此阈值，整条序列将被丢弃
- **取值范围**: 0-100
- **默认值**: 40

| 设置值 | 效果 | 适用场景 |
|--------|------|----------|
| 20-30 | 严格过滤 | 高质量要求 |
| 40-50 | 标准过滤 | 一般质控需求 |
| 60-80 | 宽松过滤 | 保留更多数据 |

### N碱基限制 (-n, --n-base-limit)

- **参数含义**: 序列中允许的最大N碱基(未知碱基)数量
- **默认值**: 10

| 设置值 | 效果 | 适用场景 |
|--------|------|----------|
| 0-5 | 严格过滤 | 对未知碱基容忍度低 |
| 10 | 标准过滤 | 一般质控需求 |
| 15-20 | 宽松过滤 | 保留更多数据 |

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **fastp** (版本 0.20.1 或更新)
  - GitHub: https://github.com/OpenGene/fastp
  - 安装方式: conda/brew/apt/yum

- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面
  - `pathlib` - 路径处理

### 安装依赖软件 | Installing Dependencies

```bash
# 使用conda安装fastp(推荐)
conda install -c bioconda fastp

# 使用Homebrew安装(macOS)
brew install fastp

# 使用apt安装(Ubuntu/Debian)
sudo apt-get install fastp

# 使用yum安装(CentOS/RHEL)
sudo yum install fastp
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器(推荐4核以上)
- **RAM**: 最少4GB(建议8GB以上)
- **存储**: 预留输入数据大小2倍的磁盘空间

## 注意事项 | Important Notes

1. **文件配对**: 双末端模式下，程序会自动查找配对文件
2. **输出覆盖**: 默认跳过已处理的样本，使用`--force`强制覆盖
3. **线程设置**: 建议根据CPU核心数设置合适的线程数
4. **内存使用**: 每个线程大约使用500MB-1GB内存
5. **文件权限**: 确保对输出目录有写入权限

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "fastp: command not found" 错误**
```bash
# 检查fastp是否安装
which fastp

# 如果未找到，使用conda安装
conda install -c bioconda fastp

# 或指定fastp路径
biopytools fastp -i input -o output --fastp-path /path/to/fastp
```

**Q: 找不到配对文件**
```bash
# 检查文件命名是否匹配
ls raw_data/*_1.fq.gz
ls raw_data/*_2.fq.gz

# 使用正确的文件后缀
biopytools fastp -i input -o output \
    --read1-suffix .R1.fq.gz \
    --read2-suffix .R2.fq.gz
```

**Q: 权限错误**
```bash
# 确保输出目录有写入权限
chmod 755 /path/to/output_dir

# 或使用有权限的目录
biopytools fastp -i input -o ~/output_data/
```

**Q: 内存不足**
```bash
# 减少线程数
biopytools fastp -i input -o output -t 4
```

**Q: 质控后数据量太少**
```bash
# 放宽质控参数
biopytools fastp -i input -o output \
    -q 20 \
    -l 30 \
    -u 50
```

**Q: 查看详细日志进行调试**
```bash
# 启用DEBUG模式
biopytools fastp -i input -o output -vv

# 保存日志到文件
biopytools fastp -i input -o output --log-file fastp.log
```

## 最佳实践 | Best Practices

### 1. 质控前检查 | Pre-Quality Control Check

```bash
# 检查原始数据
ls -lh raw_data/
find raw_data -name "*_1.fq.gz" | wc -l

# 使用fastqc查看质量
fastqc raw_data/*.fq.gz -o qc_reports/
```

### 2. 参数选择建议 | Parameter Selection Guide

| 应用场景 | 质量阈值 | 最小长度 | 不合格碱基% | N碱基限制 |
|----------|----------|----------|-------------|-----------|
| RNA-seq | 25-30 | 50-75 | 40-50 | 10 |
| WGS | 30-35 | 50-100 | 30-40 | 10 |
| 16S rRNA | 25 | 200+ | 50 | 5 |
| ChIP-seq | 20-25 | 30-50 | 50 | 15 |

### 3. 质控后验证 | Post-Quality Control Validation

```bash
# 检查清洁数据
ls -lh clean_data/*.clean.fq.gz

# 再次使用fastqc验证质控效果
fastqc clean_data/*.clean.fq.gz -o qc_reports_after/

# 使用multiqc汇总质控报告
multiqc clean_data/fastp_reports/ -o multiqc_report/
```

### 4. 性能优化 | Performance Optimization

```bash
# 根据CPU核心数设置线程
biopytools fastp -i input -o output -t $(nproc)

# 设置为核心数的80%
biopytools fastp -i input -o output -t $(($(nproc) * 4 / 5))
```

## 结果解读指南 | Result Interpretation Guide

### 质控报告关键指标 | Key Quality Control Metrics

| 指标 | 含义 | 优秀标准 |
|------|------|----------|
| **平均质量** | 所有碱基的平均质量分数 | >30 |
| **Q30比例** | 质量值>=30的碱基比例 | >85% |
| **GC含量** | G和C碱基的百分比 | 根据物种而定 |
| **接头含量** | 检测到的接头序列比例 | <5% |
| **过滤率** | 被过滤掉的reads比例 | <20% |

### HTML报告查看指南 | HTML Report Viewing Guide

1. **打开HTML报告**: 在浏览器中打开`fastp_reports/*.html`
2. **质量分布**: 查看质量随位置的变化趋势
3. **GC含量**: 与预期GC含量对比
4. **长度分布**: 检查过滤前后的长度变化
5. **接头检测**: 确认接头已被有效去除

## 相关资源 | Related Resources

- [fastp官方文档](https://github.com/OpenGene/fastp)
- [fastp使用教程](https://github.com/OpenGene/fastp/wiki)
- [FASTQ格式说明](https://en.wikipedia.org/wiki/FASTQ_format)
- [Phred质量分数](https://en.wikipedia.org/wiki/Phred_quality_score)

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用fastp相关文献：

```
Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018).
fastp: an ultra-fast all-in-one FASTQ preprocessor.
Bioinformatics, 34(17), i884-i890.
```
