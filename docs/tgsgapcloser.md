# TGS-GapCloser Gap填充分析模块

**使用三代测序数据填充基因组组装中的Gap | Fill Gaps in Genome Assembly Using TGS Long Reads**

## 功能概述 | Overview

TGS-GapCloser Gap填充分析模块是一个使用三代测序数据（ONT/PacBio）填充基因组组装中Gap的工具。该模块提供了完整的Gap填充流程，支持多种纠错模式（Racon、Pilon），自动参数设置，以及详细的日志记录，适用于各种基因组组装的质量提升工作。

## 主要特性 | Key Features

- **多平台支持**: 支持ONT、PacBio、HiFi三种三代测序数据类型
- **自动参数设置**: 根据数据类型自动设置最小同一性和最小匹配长度参数
- **多种纠错模式**: 支持none（预校正）、racon（长读长校正）、pilon（短读长校正）三种纠错模式
- **高效Gap填充**: 使用分块处理策略，支持大基因组组装
- **灵活过滤参数**: 可配置的过滤参数（reads数量、候选数等）优化填充效果
- **详细日志记录**: 完整的处理过程日志，便于问题追踪和调试
- **质量控制**: 内置输出文件验证，确保结果可靠性

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# ONT数据，无纠错模式（预校正reads）
biopytools gap-fill \
    -s scaffolds.fa \
    -t ont \
    -ir ont_reads.fa \
    -o output

# PacBio数据，使用Racon纠错
biopytools gap-fill \
    -s scaffolds.fa \
    -t pb \
    -ir pacbio_reads.fa \
    -m racon \
    -racon /path/to/racon \
    -o output

# HiFi数据，使用Pilon纠错
biopytools gap-fill \
    -s scaffolds.fa \
    -t hifi \
    -ir hifi_reads.fa \
    -m pilon \
    -pilon /path/to/pilon \
    -ngs illumina_reads.fa \
    -samtools /path/to/samtools \
    -java /path/to/java \
    -o output
```

### 高级用法 | Advanced Usage

```bash
# 自定义过滤参数和线程数
biopytools gap-fill \
    -s scaffolds.fa \
    -t ont \
    -ir ont_reads.fa \
    -idy 0.4 \
    -l 500 \
    -threads 32 \
    -chunk 5 \
    -o output

# 使用Gap大小差异检查
biopytools gap-fill \
    -s scaffolds.fa \
    -t ont \
    -ir ont_reads.fa \
    -g-check \
    -o output

# 自定义minimap2参数
biopytools gap-fill \
    -s scaffolds.fa \
    -t ont \
    -ir ont_reads.fa \
    -minmap-arg "-x ava-ont" \
    -o output
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-s, --scaff-file` | 输入scaffold文件路径（FASTA格式） | `-s scaffolds.fa` |
| `-t, --tgstype` | TGS数据类型：ont/pb/hifi | `-t ont` |
| `-ir, --reads-file` | 输入TGS reads文件路径（FASTQ格式） | `-ir reads.fa` |
| `-o, --output-prefix` | 输出前缀 | `-o output` |

### 纠错模式配置 | Error Correction Mode Configuration

| 参数 | 描述 | 默认值 | 示例 |
|------|------|--------|------|
| `-m, --mode` | 纠错模式：none/racon/pilon | `none` | `-m racon` |
| `--tgsgapcloser-path` | TGS-GapCloser可执行文件路径 | 见下方 | `--tgsgapcloser-path /path/to/tgsgapcloser` |

**默认TGS-GapCloser路径**: `/share/org/YZWL/yzwl_lixg/software/TGS-GapCloser/tgsgapcloser`

### 过滤参数 | Filter Parameters

| 参数 | 描述 | 默认值 | 示例 |
|------|------|--------|------|
| `-idy, --min-idy` | 最小同一性 | 自动设置 | `-idy 0.3` |
| `-l, --min-match` | 最小匹配长度 | 自动设置 | `-l 300` |
| `-min-nread` | 最小reads数量 | `1` | `-min-nread 5` |
| `-max-nread` | 最大reads数量（-1表示无限制） | `-1` | `-max-nread 100` |
| `-max-candidate` | 最大候选数 | `200` | `-max-candidate 500` |

**自动参数设置规则**:
- ONT: min_idy=0.3, min_match=300
- PB/HiFi: min_idy=0.2, min_match=200

### 处理配置 | Processing Configuration

| 参数 | 描述 | 默认值 | 示例 |
|------|------|--------|------|
| `-threads, --threads` | 线程数 | `16` | `-threads 32` |
| `-chunk` | 分块数量 | `3` | `-chunk 5` |
| `-g-check` | 启用Gap大小差异检查 | `False` | `-g-check` |
| `-minmap-arg` | 自定义minimap2参数 | 无 | `-minmap-arg "-x ava-ont"` |

### Racon纠错参数 | Racon Error Correction Parameters

| 参数 | 描述 | 默认值 | 示例 |
|------|------|--------|------|
| `-racon, --racon-path` | Racon可执行文件路径 | 无 | `-racon /path/to/racon` |
| `-racon-round` | Racon迭代轮数 | `3` | `-racon-round 5` |

### Pilon纠错参数 | Pilon Error Correction Parameters

| 参数 | 描述 | 默认值 | 示例 |
|------|------|--------|------|
| `-pilon, --pilon-path` | Pilon JAR文件路径 | 无 | `-pilon /path/to/pilon.jar` |
| `-ngs, --ngs-file` | NGS reads文件（BAM格式） | 无 | `-ngs illumina.bam` |
| `-java, --java-path` | Java可执行文件路径 | 无 | `-java /path/to/java` |
| `-samtools, --samtools-path` | Samtools可执行文件路径 | 无 | `-samtools /path/to/samtools` |
| `-pilon-mem` | Pilon内存分配 | `300G` | `-pilon-mem 500G` |
| `-pilon-round` | Pilon迭代轮数 | `3` | `-pilon-round 5` |

## 输出文件 | Output Files

### 主要输出 | Primary Output

| 文件 | 描述 |
|------|------|
| `{output_prefix}.gapcloser.fa` | Gap填充后的基因组序列 |
| `tgsgapcloser.log` | 详细的运行日志 |

### 纠错模式输出 | Error Correction Mode Output

使用racon或pilon模式时，还会生成纠错后的中间文件。

## 使用场景 | Use Cases

### 场景1：ONT数据Gap填充

```bash
# 使用ONT数据填充scaffold中的Gap
biopytools gap-fill \
    -s ec assemblies.fasta \
    -t ont \
    -ir ont_reads.fasta \
    -threads 24 \
    -o ec_gapfilled
```

### 场景2：PacBio数据配合Racon纠错

```bash
# PacBio数据，使用Racon进行纠错
biopytools gap-fill \
    -s assembly.fasta \
    -t pb \
    -ir pacbio_reads.fasta \
    -m racon \
    -racon /opt/racon/racon \
    -racon-round 5 \
    -threads 32 \
    -o assembly_racon_gapfilled
```

### 场景3：HiFi数据配合Pilon纠错

```bash
# HiFi数据，使用Pilon和NGS数据进行联合纠错
biopytools gap-fill \
    -s assembly.fasta \
    -t hifi \
    -ir hifi_reads.fasta \
    -m pilon \
    -pilon /opt/pilon/pilon.jar \
    -ngs illumina.bam \
    -java /opt/java/java \
    -samtools /opt/samtools/samtools \
    -pilon-mem 500G \
    -threads 32 \
    -o assembly_pilon_gapfilled
```

### 场景4：严格过滤参数

```bash
# 使用更严格的过滤参数提高填充质量
biopytools gap-fill \
    -s assembly.fasta \
    -t ont \
    -ir ont_reads.fasta \
    -idy 0.4 \
    -l 500 \
    -min-nread 10 \
    -max-candidate 100 \
    -g-check \
    -o assembly_strict
```

## 技术细节 | Technical Details

### Gap填充原理 | Gap Filling Principle

TGS-GapCloser通过以下步骤填充Gap：

1. **识别Gap位点**: 在scaffold中识别N-gap区域
2. **reads比对**: 使用minimap2将TGS reads比对到Gap侧翼序列
3. **候选reads选择**: 根据同一性和匹配长度筛选高质量reads
4. **Gap序列构建**: 使用选定的reads构建Gap填充序列
5. **纠错处理**（可选）: 使用Racon或Pilon进行纠错
6. **结果输出**: 输出填充后的基因组序列

### 纠错模式选择 | Error Correction Mode Selection

| 模式 | 适用场景 | 优点 | 缺点 |
|------|----------|------|------|
| **none** | 已校正的高质量reads | 速度快，资源消耗低 | 填充准确性依赖输入数据质量 |
| **racon** | 原始长读长数据 | 提高长读长准确性 | 需要额外计算资源 |
| **pilon** | 有NGS数据可用 | 短读长精度高，纠错效果好 | 需要NGS数据和更多资源 |

### 参数优化建议 | Parameter Optimization Recommendations

**ONT数据**:
- min_idy: 0.3-0.4（默认0.3）
- min_match: 300-500（默认300）
- 推荐纠错模式: racon

**PacBio数据**:
- min_idy: 0.2-0.3（默认0.2）
- min_match: 200-300（默认200）
- 推荐纠错模式: racon

**HiFi数据**:
- min_idy: 0.2-0.3（默认0.2）
- min_match: 200-300（默认200）
- 推荐纠错模式: none（HiFi数据已较准确）

## 常见问题 | FAQ

**Q: 如何选择合适的纠错模式？**

A:
- 如果输入数据已经过校正（如CCS、HiFi），使用`none`模式
- 如果是原始ONT/PacBio数据，推荐使用`racon`模式
- 如果有高质量的NGS数据，可以使用`pilon`模式获得最佳效果

**Q: 内存不足怎么办？**

A:
- 减少线程数: `-threads 8`
- 增加分块数: `-chunk 5`
- 对于Pilon模式，减少内存分配: `-pilon-mem 100G`

**Q: Gap填充效果不佳？**

A:
- 检查reads质量，确保数据覆盖度足够
- 调整过滤参数，尝试更严格的设置（提高min_idy和min_match）
- 启用Gap大小差异检查: `-g-check`
- 考虑使用纠错模式提高准确性

**Q: 如何加速处理？**

A:
- 增加线程数: `-threads 32`
- 减少分块数: `-chunk 2`
- 使用预校正reads，选择`none`模式

## 依赖要求 | Dependencies

### 必需依赖 | Required Dependencies

- TGS-GapCloser: 默认路径 `/share/org/YZWL/yzwl_lixg/software/TGS-GapCloser/tgsgapcloser`
- Python 3.7+

### 可选依赖 | Optional Dependencies

- Racon: 用于`racon`纠错模式
- Pilon: 用于`pilon`纠错模式
- Java: Pilon运行需要
- Samtools: Pilon模式处理BAM文件需要

## 引用 | Citation

如果使用TGS-GapCloser模块，请引用：

> TGS-GapCloser: Filling Gaps in Genome Assembly Using TGS Long Reads
> DOI: [待添加]

## 版本历史 | Version History

- **v1.0.0** (2025): 初始版本，支持基本Gap填充功能
  - 支持ONT、PacBio、HiFi数据类型
  - 支持none、racon、pilon三种纠错模式
  - 自动参数设置和分块处理
  - 详细的日志记录和输出验证
