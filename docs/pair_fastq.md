# 📝 FASTQ配对修复分析模块

**专业的配对混乱FASTQ文件修复工具 | Professional Tool for Fixing Paired-end FASTQ Files with Pairing Issues**

## 📖 功能概述 | Overview

FASTQ配对修复分析模块是一个专门用于修复配对混乱的双端测序FASTQ文件的工具，支持**seqkit**和**repair.sh (BBMap)**两种工具，提供自动化的配对检测、修复和验证功能。支持批量处理、灵活的文件命名规则配置和详细的日志记录，适用于各种高通量测序数据的配对问题修复。

## ✨ 主要特性 | Key Features

- **🔍 智能配对检测**: 自动识别并匹配R1和R2文件
- **🔄 双工具支持**: 支持seqkit和repair.sh (BBMap)两种处理引擎
- **📂 智能文件组织**: 自动创建paired和single子目录，分离配对和单条reads
- **🔄 批量处理**: 一次性处理多个样本的配对修复
- **⚙️ 灵活命名规则**: 支持自定义文件名后缀模式
- **🛡️ 完整性验证**: 自动检测缺失的R1或R2文件
- **📊 详细统计**: 提供完整的处理结果和错误报告
- **🚀 高效处理**: 支持多线程并行处理
- **📝 完善日志**: 记录所有处理步骤和错误信息

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本配对修复（默认使用repair.sh）
biopytools pair-fastq \
    -i /path/to/raw_data \
    -o /path/to/fixed_data

# 使用多线程加速
biopytools pair-fastq \
    -i /path/to/raw_data \
    -o /path/to/fixed_data \
    -t 16

# 使用seqkit工具（适合小文件）
biopytools pair-fastq \
    -i /path/to/raw_data \
    -o /path/to/fixed_data \
    --tool seqkit
```

### 高级用法 | Advanced Usage

```bash
# 自定义文件后缀和输出路径
biopytools pair-fastq \
    -i raw_sequences \
    -o corrected_sequences \
    --suffix1 ".R1.fastq.gz" \
    --suffix2 ".R2.fastq.gz" \
    -t 32 \
    --verbose

# 使用repair.sh并指定内存
biopytools pair-fastq \
    -i raw_data \
    -o fixed_data \
    --tool repair \
    --repair-memory 500g \
    -t 32

# Dry run模式（仅显示命令）
biopytools pair-fastq \
    -i raw_data \
    -o fixed_data \
    --dry-run
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入目录（包含FASTQ文件）| `-i ./raw_data` |
| `-o, --output` | 输出目录 | `-o ./fixed_data` |

### 处理参数 | Processing Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `12` | 🔧 线程数（用于并行处理）|
| `--suffix1` | `_1.fq.gz` | 📄 R1文件后缀模式 |
| `--suffix2` | `_2.fq.gz` | 📄 R2文件后缀模式 |
| `--tool` | `repair` | 🛠️ 工具选择（seqkit或repair） |

### 工具配置 | Tool Configuration

**seqkit参数**:

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--seqkit-bin` | `seqkit` | 🛠️ seqkit二进制文件路径 |

**repair.sh参数**:

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--repair-sh` | `repair.sh` | 🔧 repair.sh脚本名称 |
| `--repair-conda-env` | `bbmap_v.39.81` | 🐍 conda环境名称 |
| `--repair-memory` | `300g` | 💾 内存参数（-Xmx）|

### 调试和日志选项 | Debug and Logging Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--dry-run` | `False` | 🎯 仅显示命令不执行 |
| `--verbose` | `False` | 📢 详细输出模式 |
| `--log-file` | 自动生成 | 📝 日志文件路径 |
| `--log-level` | `INFO` | 📊 日志级别（DEBUG/INFO/WARNING/ERROR/CRITICAL）|

## 📁 输入文件格式 | Input File Formats

### 文件命名规则 | File Naming Rules

工具通过文件名后缀自动识别配对：

```
示例1：默认后缀
sampleA_1.fq.gz  ← R1文件
sampleA_2.fq.gz  ← R2文件

示例2：自定义后缀
sampleA.R1.fastq.gz  ← R1文件
sampleA.R2.fastq.gz  ← R2文件
```

### 输入目录结构 | Input Directory Structure

```
raw_data/
├── sample1_1.fq.gz
├── sample1_2.fq.gz
├── sample2_1.fq.gz
├── sample2_2.fq.gz
├── sample3_1.fq.gz
└── sample3_2.fq.gz
```

## 💡 使用示例 | Usage Examples

### 示例1：基本配对修复 | Example 1: Basic Pair Fixing

```bash
# 修复配对混乱的FASTQ文件
biopytools pair-fastq \
    -i ./raw_sequences \
    -o ./fixed_sequences \
    -t 16
```

**适用场景**：
- 测序数据由于某些原因导致配对混乱
- 文件名正确但内部reads顺序不匹配
- 需要快速修复大批量样本

### 示例2：自定义文件后缀 | Example 2: Custom File Suffix

```bash
# 处理不同命名规则的文件
biopytools pair-fastq \
    -i ./data \
    -o ./corrected \
    --suffix1 ".R1.fastq.gz" \
    --suffix2 ".R2.fastq.gz" \
    -t 24
```

**适用场景**：
- 测序公司使用非标准命名规则
- 使用了特殊的文件名格式
- 需要处理来自不同来源的数据

### 示例3：Dry Run测试 | Example 3: Dry Run Testing

```bash
# 查看将要执行的命令但不实际运行
biopytools pair-fastq \
    -i ./test_data \
    -o ./output \
    --dry-run \
    --verbose
```

**适用场景**：
- 测试和验证处理流程
- 检查命令是否正确
- 估算处理时间

### 示例4：详细日志记录 | Example 4: Detailed Logging

```bash
# 启用详细日志和自定义日志文件
biopytools pair-fastq \
    -i ./large_dataset \
    -o ./fixed_dataset \
    -t 32 \
    --verbose \
    --log-file ./pair_fixing.log \
    --log-level DEBUG
```

**适用场景**：
- 处理大规模数据集
- 需要详细的调试信息
- 问题排查和性能分析

### 示例5：处理特定样本子集 | Example 5: Process Specific Subset

```bash
# 通过子目录处理特定样本
biopytools pair-fastq \
    -i ./raw_data/samples_subset \
    -o ./fixed_data/samples_subset \
    -t 16 \
    --log-level WARNING
```

**适用场景**：
- 只处理部分样本
- 分批处理大规模数据
- 测试和验证阶段

## 📊 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

工具会自动创建`paired`和`single`两个子目录来组织输出文件：

```
fixed_data/
├── paired/                              # 配对的R1和R2文件
│   ├── sample1_1.fq.gz                  # 修复后的R1文件
│   ├── sample1_2.fq.gz                  # 修复后的R2文件
│   ├── sample2_1.fq.gz
│   ├── sample2_2.fq.gz
│   ├── sample3_1.fq.gz
│   └── sample3_2.fq.gz
├── single/                              # 单条reads文件
│   ├── sample1_single.fq.gz             # 未配对的reads
│   ├── sample2_single.fq.gz
│   └── sample3_single.fq.gz
└── logs/
    └── pair_fastq_20260404_120000.log    # 处理日志
```

**目录说明**：
- **paired/**: 存放成功配对的R1和R2文件，文件名与输入保持一致
- **single/**: 存放未配对的单条reads（repair.sh必需输出，seqkit自动生成）
- **logs/**: 存放处理日志文件

### 日志文件内容 | Log File Content

日志文件包含：
- ✅ 配置信息记录
- ✅ 样本发现和配对统计
- ✅ 每个样本的处理状态
- ✅ 错误和警告信息
- ✅ 最终处理结果汇总

### 控制台输出示例 | Console Output Example

```
============================================================
开始FASTQ配对修复流程|Starting FASTQ pair fixing pipeline
============================================================
配置信息|Configuration:
  输入目录|Input directory: /path/to/raw_data
  输出目录|Output directory: /path/to/fixed_data
  配对文件|Paired files: /path/to/fixed_data/paired/
  单条reads|Singleton reads: /path/to/fixed_data/single/
  R1后缀|R1 suffix: _1.fq.gz
  R2后缀|R2 suffix: _2.fq.gz
  工具|Tool: repair
  线程数|Threads: 16
  Repair内存|Repair memory: 300g
------------------------------------------------------------
扫描输入目录: /path/to/raw_data|Scanning input directory
找到 50 对完整文件|Found 50 complete pairs
开始处理 50 个样本|Start processing 50 samples
------------------------------------------------------------
[1/50] 处理中...|Processing...
处理样本: sample1|Processing sample: sample1
  样本 sample1 完成|Sample sample1 completed
  单条reads|Singleton reads: /path/to/fixed_data/single/sample1_single.fq.gz (12.34 MB)
...
============================================================
处理结果汇总|Processing Summary:
  总样本数|Total samples: 50
  成功|Success: 50/50
  失败|Failed: 0/50
============================================================
```

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

**处理工具（二选一）**:

- **seqkit** (版本 2.0+) - 适合小文件
  - 下载地址: https://github.com/shenwei356/seqkit
  - 开源工具，无需许可证
  - 轻量级，速度快

- **repair.sh (BBMap)** (版本 38+) - 适合大文件（默认）
  - 包含在BBMap套件中
  - 需要conda环境支持
  - 内存需求大，处理能力更强

**Python环境**:

- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面
  - `pathlib` - 路径处理
  - `dataclasses` - 配置管理
  - `subprocess` - 系统调用

### 安装依赖软件 | Installing Dependencies

```bash
# 安装seqkit（可选）
# 方法1：使用conda（推荐）
conda install -c bioconda seqkit

# 方法2：从源码编译
go install github.com/shenwei356/seqkit/v2/cmd/seqkit@latest

# 方法3：下载预编译二进制文件
wget https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz
tar -xzf seqkit_linux_amd64.tar.gz
chmod +x seqkit
sudo mv seqkit /usr/local/bin/

# 验证seqkit安装
seqkit version

# 安装repair.sh (BBMap)
# 使用conda安装BBMap
conda create -n bbmap_v.39.81 -c bioconda bbmap=39.81
conda activate bbmap_v.39.81

# 验证repair.sh安装
conda run -n bbmap_v.39.81 repair.sh --help
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐12核以上以充分利用repair.sh）
- **RAM**:
  - 使用seqkit: 最少4GB（大规模数据集推荐16GB以上）
  - 使用repair.sh: 推荐至少300GB（可调整`--repair-memory`参数）
- **存储**: 预留与输入文件大小相同的磁盘空间
- **I/O**: SSD硬盘可显著提升处理速度

## ⚠️ 注意事项 | Important Notes

### 1. 文件配对要求 | File Pairing Requirements

- ✅ R1和R2文件必须有相同的样本名称（前缀）
- ✅ 文件名必须符合指定的后缀模式
- ✅ 每个样本只能有一个R1和一个R2文件

### 2. 工具选择建议 | Tool Selection Recommendations

- **seqkit**: 适合小文件（<10GB）、快速测试、资源受限环境
- **repair.sh**: 适合大文件（>10GB）、大规模生产环境、需要高稳定性

### 3. 处理前检查 | Pre-processing Checks

```bash
# 检查输入目录
ls -lh /path/to/raw_data

# 验证seqkit可用
seqkit version

# 验证repair.sh可用
conda run -n bbmap_v.39.81 repair.sh --help

# 查看磁盘空间
df -h /path/to/output
```

### 4. 处理建议 | Processing Recommendations

- ⚠️ 建议先用少量样本测试参数
- ⚠️ 大规模数据推荐使用repair.sh和更多线程
- ⚠️ 保留原始数据备份
- ⚠️ 使用日志文件记录处理过程
- ⚠️ 注意paired和single子目录的文件组织结构

### 5. 常见问题 | Common Issues

❌ **问题1：找不到配对文件**
```
错误|Error: 找到 0 对完整文件|Found 0 complete pairs
```
**解决方案**：
- 检查文件命名是否符合后缀规则
- 使用`--suffix1`和`--suffix2`调整后缀模式

❌ **问题2：工具未找到**
```
错误|Error: seqkit工具未找到或不可执行
错误|Error: conda环境未找到或不可用
```
**解决方案**：
- seqkit: `conda install -c bioconda seqkit`
- repair.sh: 确认conda环境存在，或使用`--repair-conda-env`指定

❌ **问题3：内存不足（repair.sh）**
```
错误|Error: Cannot allocate memory
```
**解决方案**：
- 减少内存分配：`--repair-memory 100g`
- 或切换到seqkit：`--tool seqkit`

## 🐛 故障排除 | Troubleshooting

### 调试模式 | Debug Mode

```bash
# 启用详细日志和调试输出
biopytools pair-fastq \
    -i ./test_data \
    -o ./output \
    --verbose \
    --log-level DEBUG \
    --log-file debug.log
```

### 常见错误和解决方案 | Common Errors and Solutions

| 错误信息 | 可能原因 | 解决方案 |
|---------|---------|---------|
| `输入目录不存在` | 路径错误 | 检查输入路径是否正确 |
| `找到 0 对完整文件` | 文件命名不符合规则 | 调整`--suffix1`和`--suffix2` |
| `seqkit工具未找到` | seqkit未安装 | `conda install -c bioconda seqkit` |
| `conda环境未找到` | conda环境不存在 | 检查`--repair-conda-env`或创建环境 |
| `内存不足` | repair.sh内存不足 | 调整`--repair-memory`或切换到seqkit |

## 📚 相关资源 | Related Resources

- [seqkit官方文档](https://bioinf.shenwei.me/seqkit/)
- [seqkit GitHub仓库](https://github.com/shenwei356/seqkit)
- [BBMap官方文档](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)
- [BBMap SourceForge](https://sourceforge.net/projects/bbmap/)
- [FASTQ格式规范](https://en.wikipedia.org/wiki/FASTQ_format)
- [双端测序原理](https://www.illumina.com/science/technology/next-generation-sequencing.html)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

---

## 🔬 引用信息 | Citation

如果在研究中使用此工具，请引用：

```bibtex
@software{biopytools_pair_fastq,
  title = {BioPyTools FASTQ Pair Fixing Module},
  author = {LI, Xiang},
  year = {2026},
  url = {https://github.com/yourusername/biopytools}
}
```

**seqkit引用**：
```
Shen, W., Le, S., Li, Y., & Hu, F. (2016).
SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation.
PLoS ONE, 11(10), e0163962.
```
