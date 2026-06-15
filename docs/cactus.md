# Cactus泛基因组分析工具 | Cactus Pangenome Analysis Tool

版本 | Version: 1.0.0
作者 | Author: Xiang LI
日期 | Date: 2026-04-02

## 概述 | Overview

Cactus工具是基于**Cactus**和**Minigraph**的泛基因组构建Python封装，使用Singularity容器调用Cactus，提供完整的泛基因组构建流程。

The Cactus tool is a Python wrapper for **Cactus** and **Minigraph** pangenome construction, using Singularity containers to invoke Cactus, providing a complete pangenome construction pipeline.

## 功能特点 | Features

- 🧬 **完整流程**: Minigraph-Cactus 5步骤流程一站式完成
- 📊 **多格式输出**: 支持GFA、HAL、GBZ、ODGI等多种输出格式
- 🎯 **Singularity容器**: 使用Singularity确保环境一致性
- 💾 **断点续传**: 支持中断恢复，自动跳过已完成步骤
- ⚡ **高性能**: 支持多线程并行计算，可配置内存限制
- 📈 **自动验证**: 自动检查Singularity、Cactus镜像和输入文件
- 🔧 **灵活配置**: 丰富的参数配置选项，适应不同数据类型

## 分析流程 | Analysis Pipeline

Cactus泛基因组构建包含5个主要步骤：

Cactus pangenome construction consists of 5 main steps:

### 步骤1: cactus-minigraph | Minigraph图构建

使用Minigraph构建初始泛基因组图。

Uses Minigraph to build initial pangenome graph.

### 步骤2: cactus-graphmap | 序列比对

将所有基因组序列比对到泛基因组图。

Aligns all genome sequences to the pangenome graph.

### 步骤3: cactus-graphmap-split | 分割比对结果

分割比对结果以进行并行处理。

Splits alignment results for parallel processing.

### 步骤4: cactus-align | 多序列比对

执行多序列比对以生成最终泛基因组。

Performs multiple sequence alignment to generate final pangenome.

### 步骤5: cactus-graphmap-join | 合并结果

合并所有比对结果并生成输出文件。

Merges all alignment results and generates output files.

## 安装和使用 | Installation and Usage

### 前置要求 | Prerequisites

#### Singularity安装 | Singularity Installation

```bash
# 使用conda安装Singularity|Install Singularity using conda
conda create -n singularity_v.3.8.7 -c conda-forge singularity
conda activate singularity_v.3.8.7

# 或从系统包管理器安装|Or install from system package manager
# Ubuntu/Debian
sudo apt-get install singularity-container

# CentOS/RHEL
sudo yum install singularity
```

#### Cactus SIF镜像 | Cactus SIF Image

```bash
# 下载Cactus SIF镜像|Download Cactus SIF image
# 方式1: 从Docker Hub转换|Method 1: Convert from Docker Hub
singularity pull cactus_v3.1.4.sif docker://quay.io/comparative-genomics-toolkit/cactus:v3.1.4

# 方式2: 从现有位置复制|Method 2: Copy from existing location
cp /path/to/existing/cactus_v3.1.4.sif /share/org/YZWL/yzwl_lixg/software/singularity/
```

#### Python环境要求 | Python Environment Requirements

```bash
# biopytools会自动处理Python依赖
# 确保biopytools已安装
pip install biopytools
```

### 基本用法 | Basic Usage

#### 准备序列文件 | Prepare Sequence File

创建序列文件（seqfile.txt），每行一个基因组文件路径：

Create sequence file (seqfile.txt), one genome path per line:

```
/path/to/reference.fa
/path/to/sample1.fa
/path/to/sample2.fa
/path/to/sample3.fa
```

**重要提示|Important Notes:**
- 第一个文件必须是参考基因组|First file must be reference genome
- 支持相对路径和绝对路径|Supports relative and absolute paths
- 建议使用绝对路径|Absolute paths recommended

#### 运行分析 | Run Analysis

```bash
# 基本用法|Basic usage
biopytools cactus -s seqfile.txt -o output/ -r reference

# 指定输出格式|Specify output formats
biopytools cactus -s seqfile.txt -o output/ -r reference --formats gfa hal odgi

# 自定义Singularity路径|Custom Singularity path
biopytools cactus -s seqfile.txt -o output/ -r reference \
  --singularity ~/miniforge3/envs/singularity_v.3.8.7/bin/singularity \
  --cactus-sif ~/images/cactus_v3.1.4.sif

# 增加线程和内存|Increase threads and memory
biopytools cactus -s seqfile.txt -o output/ -r reference \
  --threads 24 --max-memory 200G

# 保留jobstore用于调试|Keep jobstore for debugging
biopytools cactus -s seqfile.txt -o output/ -r reference --no-cleanup
```

### 参数说明 | Parameter Description

#### 必需参数 | Required Arguments

| 参数 | 说明 | 示例 |
|------|------|------|
| `-s, --seqfile` | 序列文件路径|Sequence file path | `-s seqfile.txt` |
| `-o, --output` | 输出目录|Output directory | `-o output/` |
| `-r, --reference` | 参考基因组名称|Reference genome name | `-r reference` |

#### 可选参数 | Optional Arguments

| 参数 | 默认值 | 说明 | 示例 |
|------|--------|------|------|
| `--singularity` | `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity` | Singularity可执行文件路径|Singularity executable path | `--singularity /usr/bin/singularity` |
| `--cactus-sif` | `/share/org/YZWL/yzwl_lixg/software/singularity/cactus_v3.1.4.sif` | Cactus SIF镜像路径|Cactus SIF image path | `--cactus-sif ~/images/cactus.sif` |
| `--jobstore` | `cactus-jobstore` | Toil jobstore目录名称|Toil jobstore directory name | `--jobstore my-jobstore` |
| `--out-name` | `cactus_output` | 输出文件前缀|Output file prefix | `--out-name my_pangenome` |
| `--no-cleanup` | `False` | 保留jobstore不删除|Keep jobstore without deleting | `--no-cleanup` |
| `--formats` | `gfa hal gbz odgi` | 输出格式|Output formats | `--formats gfa hal` |
| `-t, --threads` | `12` | CPU核心数|Number of CPU cores | `--threads 24` |
| `-m, --max-memory` | `100G` | 最大内存|Maximum memory | `--max-memory 200G` |
| `--bind` | `自动检测|Auto-detect` | 绑定目录到容器|Bind directory to container | `--bind /data:/data` |
| `--log-level` | `INFO` | 日志级别|Log level | `--log-level DEBUG` |

## 输出文件 | Output Files

### 主要输出文件 | Main Output Files

根据指定的输出格式，生成相应的文件：

Based on specified output formats, generates corresponding files:

- `{out_name}.gfa` - GFA格式泛基因组图|GFA format pangenome graph
- `{out_name}.hal` - HAL格式多序列比对|HAL format multiple sequence alignment
- `{out_name}.gbz` - GBZ格式压缩图|GBZ format compressed graph
- `{out_name}.odgi` - ODGI格式图|ODGI format graph

### 日志文件 | Log Files

- `{out_name}.log` - 完整运行日志|Complete run log

### 临时文件 | Temporary Files

- `cactus-jobstore/` - Toil jobstore目录（默认分析完成后删除）|Toil jobstore directory (deleted by default after completion)

## 断点续传功能 | Checkpoint Resume Feature

本工具支持断点续传，自动跳过已完成的步骤：

This tool supports checkpoint resume, automatically skipping completed steps:

### 检查机制 | Check Mechanism

检查所有指定的输出文件是否存在：

Checks if all specified output files exist:

- GFA文件（如果指定）|GFA file (if specified)
- HAL文件（如果指定）|HAL file (if specified)
- GBZ文件（如果指定）|GBZ file (if specified)
- ODGI文件（如果指定）|ODGI file (if specified)

### 使用示例 | Usage Example

```bash
# 第一次运行 - 执行完整流程
biopytools cactus -s seqfile.txt -o output/ -r reference
# 输出：开始Cactus分析|Starting Cactus analysis

# 中断后重新运行 - 自动跳过已完成步骤
biopytools cactus -s seqfile.txt -o output/ -r reference
# 输出：分析已完成，跳过|Analysis already completed, skipping
```

## 常见问题 | FAQ

### 1. Singularity路径找不到怎么办？

**问题|Question**: Singularity路径找不到怎么办？

**答|Answer**:
```bash
# 查找Singularity安装位置|Find Singularity installation location
which singularity
# 或|or
conda run -n singularity_v.3.8.7 which singularity

# 使用完整路径|Use full path
biopytools cactus -s seqfile.txt -o output/ -r reference \
  --singularity /full/path/to/singularity
```

### 2. 如何获取Cactus SIF镜像？

**问题|Question**: 如何获取Cactus SIF镜像？

**答|Answer**:
```bash
# 从Docker Hub转换|Convert from Docker Hub
singularity pull cactus_v3.1.4.sif docker://quay.io/comparative-genomics-toolkit/cactus:v3.1.4

# 或使用已有镜像|Or use existing image
cp /path/to/cactus.sif /share/org/YZWL/yzwl_lixg/software/singularity/cactus_v3.1.4.sif
```

### 3. 内存不足怎么办？

**问题|Question**: 内存不足怎么办？

**答|Answer**:
```bash
# 增加内存限制|Increase memory limit
biopytools cactus -s seqfile.txt -o output/ -r reference --max-memory 200G

# 减少线程数（降低内存使用）|Reduce threads (lower memory usage)
biopytools cactus -s seqfile.txt -o output/ -r reference --threads 8
```

### 4. 如何选择输出格式？

**问题|Question**: 如何选择输出格式？

**答|Answer**:

| 格式 | 适用场景|Use case | 推荐工具|Recommended tools |
|------|----------|-----------|-------------------|
| GFA | 可视化|Visualization | Bandage, VG | ✓ 推荐|Recommended |
| HAL | 比较分析|Comparative analysis | HAL tools | ✓ 推荐|Recommended |
| GBZ | 压缩存储|Compressed storage | ODGI | ✓ 推荐|Recommended |
| ODGI | 图分析|Graph analysis | ODGI | ✓ 推荐|Recommended |
| VG | 变异分析|Variant analysis | VG | 可选|Optional |
| PSA | 序列提取|Sequence extraction | HAL | 可选|Optional |

### 5. 参考基因组名称如何填写？

**问题|Question**: 参考基因组名称如何填写？

**答|Answer**:

参考基因组名称必须与序列文件中第一个基因组的名称匹配：

Reference genome name must match the name of the first genome in sequence file:

```bash
# seqfile.txt内容|seqfile.txt content:
/path/to/genomes/href.fa
/path/to/genomes/sample1.fa
/path/to/genomes/sample2.fa

# 运行命令|Run command:
biopytools cactus -s seqfile.txt -o output/ -r href
```

**注意|Note**: 使用文件名（不含路径和扩展名）|Use filename (without path and extension)

### 6. 如何绑定额外的目录到容器？

**问题|Question**: 如何绑定额外的目录到容器？

**答|Answer**:
```bash
# 绑定单个目录|Bind single directory
biopytools cactus -s seqfile.txt -o output/ -r reference \
  --bind /extra/data:/extra/data

# 绑定多个目录|Bind multiple directories
biopytools cactus -s seqfile.txt -o output/ -r reference \
  --bind /data:/data \
  --bind /reference:/reference
```

### 7. 分析需要多长时间？

**问题|Question**: 分析需要多长时间？

**答|Answer**:

| 基因组数量|Genome count | 基因组大小|Genome size | CPU核心数|CPU cores | 预计时间|Estimated time |
|-------------|--------------|-------------|----------|--------------|
| 3 | 1 Gb | 12 | 2-4 小时|hours |
| 10 | 1 Gb | 24 | 6-12 小时|hours |
| 50 | 1 Gb | 48 | 24-48 小时|hours |

**注意|Note**: 实际时间取决于基因组复杂度和变异程度|Actual time depends on genome complexity and variation level

## 输出文件解读 | Output File Interpretation

### GFA文件 | GFA File

GFA (Graphical Fragment Assembly) 格式描述泛基因组图结构：

GFA (Graphical Fragment Assembly) format describes pangenome graph structure:

```bash
# 使用Bandage可视化|Visualize using Bandage
Bandage image cactus_output.gfa output.png

# 使用VG工具分析|Analyze using VG tools
vg view cactus_output.gfa > cactus_output.vg
vg stats cactus_output.vg
```

### HAL文件 | HAL File

HAL (Hierarchical Alignment) 格式用于多序列比对：

HAL (Hierarchical Alignment) format for multiple sequence alignment:

```bash
# 查看HAL统计信息|View HAL statistics
halStats cactus_output.hal

# 提取特定基因组|Extract specific genome
hal2fasta cactus_output.hal sample1 > sample1.fa
```

### ODGI文件 | ODGI File

ODGI格式用于高效的图分析：

ODGI format for efficient graph analysis:

```bash
# 查看图统计信息|View graph statistics
odgi stats -i cactus_output.odgi

# 可视化图|Visualize graph
odgi draw -i cactus_output.odgi -o cactus.png
```

## 性能优化建议 | Performance Optimization Recommendations

### 1. 内存优化 | Memory Optimization

```bash
# 根据可用内存调整|Adjust based on available memory
--max-memory 200G  # 200GB内存|200GB memory
```

### 2. 速度优化 | Speed Optimization

```bash
# 增加线程数|Increase threads
--threads 48  # 使用48个CPU核心|Use 48 CPU cores
```

### 3. 存储优化 | Storage Optimization

```bash
# 只生成需要的格式|Generate only needed formats
--formats gfa hal  # 只生成GFA和HAL|Only generate GFA and HAL
```

## 依赖项 | Dependencies

### 系统依赖 | System Dependencies

- Singularity >= 3.0
- Cactus SIF镜像|Cactus SIF image (v3.1.4推荐|recommended)

### Python依赖 | Python Dependencies

- Python >= 3.8
- pathlib
- dataclasses
- logging
- subprocess

## 参考资源 | References

### Cactus相关

- **Cactus GitHub**: https://github.com/ComparativeGenomicsToolkit/cactus
- **Cactus文档**: https://cactus-genome.github.io/
- **Minigraph论文**: https://www.nature.com/articles/s41587-020-00750-5

### 泛基因组分析方法

- **泛基因组综述**: https://www.nature.com/articles/s41579-020-00452-w
- **GFA格式规范**: https://github.com/GFA-spec/GFA-spec/
- **ODGI工具**: https://github.com/vgteam/odgi

## 更新日志 | Changelog

### v1.0.0 (2026-04-02)

- 初始版本发布
- 支持Singularity容器调用Cactus
- 支持GFA、HAL、GBZ、ODGI输出格式
- 实现断点续传功能
- 完整的日志分离（stdout/stderr）
- 自动环境检测和验证
- 结果文件自动检查

## 许可证 | License

本工具遵循biopytools项目的许可证。

This tool follows the biopytools project license.

## 联系方式 | Contact

作者|Author: Xiang LI <lixiang117423@gmail.com>

项目地址|Project: https://github.com/your-org/biopytools
