# Swave 结构变异检测分析模块

**泛基因组结构变异检测工具 | Pangenome Structural Variant Detection Tool**

## 功能概述 | Overview

Swave 结构变异检测分析模块是一个基于泛基因组图的结构变异（SV）和复杂SV检测工具，支持从多个基因组组装构建的泛基因图中精确检测结构变异。支持多种泛基因组图构建工具（minigraph、minigraph-cactus、pggb），提供完整的SV检测流程和灵活的参数配置。

## 主要特性 | Key Features

- **多种图来源支持**: 支持 minigraph、minigraph-cactus 和 pggb 构建的 GFA 文件
- **高精度检测**: 基于LSTM模型的结构变异预测
- **复杂SV检测**: 支持检测复杂的结构变异组合
- **灵活样本选择**: 支持单样本、多样本和指定样本分析
- **质量分级**: 提供三级质量标签（HighQual、MediumQual、LowQual）
- **可配置参数**: 支持SV大小范围、组件数量等多维度参数调整
- **多个子命令**: 提供完整的工具链（检测、转换、提取）

## 快速开始 | Quick Start

### 查看所有子命令 | View All Subcommands

```bash
biopytools swave -h
```

### 基本用法 | Basic Usage

```bash
# 1. SV检测 - 最常用
biopytools swave call \
    -i assemblies.tsv \
    -r reference.fa \
    -g pangenome.gfa \
    -s minigraph \
    -o swave_results

# 2. 转换图路径为序列
biopytools swave convert-seq \
    --vcf-path input.vcf \
    --gfa-path pangenome.gfa \
    --ref-path reference.fa

# 3. 提取特定样本
biopytools swave extract-sample \
    --vcf-path input.vcf \
    --spec-samples sample1 sample2

# 4. 提取CSV格式
biopytools swave extract-csv \
    --vcf-path input.vcf \
    --spec-csv INV
```

## 子命令详解 | Subcommands Details

### 1. call - 检测结构变异

这是Swave的主要功能，用于从泛基因组图检测结构变异。

#### 基本用法

```bash
biopytools swave call \
    -i assemblies.tsv \
    -r reference.fa \
    -g pangenome.gfa \
    -s minigraph \
    -o output_dir
```

#### Minigraph图

```bash
biopytools swave call \
    -i assemblies.tsv \
    -r reference.fa \
    -g pangenome.gfa \
    -s minigraph \
    -o swave_results
```

#### Cactus图（需要decomposed VCF）

```bash
biopytools swave call \
    -i assemblies.tsv \
    -r reference.fa \
    -g pangenome.cactus.gfa \
    -s cactus \
    --decomposed-vcf decomposed.vcf.gz \
    -o swave_results
```

#### 指定样本检测

```bash
biopytools swave call \
    -i assemblies.tsv \
    -r reference.fa \
    -g pangenome.gfa \
    -s minigraph \
    --spec-samples sample1 sample2 sample3 \
    -o specific_samples
```

#### 自定义SV大小和质量控制

```bash
biopytools swave call \
    -i assemblies.tsv \
    -r reference.fa \
    -g pangenome.gfa \
    -s minigraph \
    --min-sv-size 100 \
    --max-sv-size 500000 \
    --max-sv-comps 3 \
    -t 8 \
    -o strict_sv_detection
```

### 2. convert-seq - 转换图路径为序列

将泛基因组图的路径转换为VCF文件的REF和ALT列的序列。

#### 基本用法

```bash
biopytools swave convert-seq \
    --vcf-path input.vcf \
    --gfa-path pangenome.gfa \
    --ref-path reference.fa
```

#### 输出为pangenie格式

```bash
biopytools swave convert-seq \
    --vcf-path input.vcf \
    --gfa-path pangenome.gfa \
    --ref-path reference.fa \
    --force-pangenie \
    --output-path pangenie_output.vcf
```

### 3. convert-plines - 转换VCF为GFA P lines

将minigraph call VCF转换为GFA P lines格式。

#### 基本用法

```bash
biopytools swave convert-plines \
    --gfa-path pangenome.gfa \
    --vcf-path input.vcf
```

#### 使用参考VCF

```bash
biopytools swave convert-plines \
    --gfa-path pangenome.gfa \
    --vcf-path input.vcf \
    --ref-vcf-path reference.vcf \
    --output-path output.gfa
```

#### 强制vg格式

```bash
biopytools swave convert-plines \
    --gfa-path pangenome.gfa \
    --vcf-path input.vcf \
    --force-vg \
    --output-path vg_format.gfa
```

### 4. extract-csv - 从VCF提取CSV

从VCF文件提取特定类型的CSV格式（INV或DUP）。

#### 提取倒位

```bash
biopytools swave extract-csv \
    --vcf-path input.vcf \
    --spec-csv INV \
    --output-path inversions.csv
```

#### 提取重复

```bash
biopytools swave extract-csv \
    --vcf-path input.vcf \
    --spec-csv DUP \
    --output-path duplications.csv
```

#### 提取所有类型

```bash
biopytools swave extract-csv \
    --vcf-path input.vcf \
    --spec-csv All \
    --output-path all_sv.csv
```

### 5. extract-sample - 提取特定样本的SV

从VCF文件中提取特定样本的结构变异。

#### 基本用法

```bash
biopytools swave extract-sample \
    --vcf-path input.vcf \
    --spec-samples sample1 sample2 \
    --output-path extracted.vcf
```

#### 单个样本

```bash
biopytools swave extract-sample \
    --vcf-path input.vcf \
    --spec-samples sample1 \
    --output-path sample1_only.vcf
```

## 参数说明 | Parameters

### call 子命令参数

#### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --assemblies-tsv` | 样本组装TSV文件 | `-i assemblies.tsv` |
| `-r, --ref-fasta` | 参考基因组FASTA文件 | `-r reference.fa` |
| `-g, --gfa-file` | 泛基因组图GFA文件 | `-g pangenome.gfa` |
| `-s, --gfa-source` | GFA来源 (minigraph/cactus/pggb) | `-s minigraph` |

#### 软件配置 | Software Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--swave-path` | `~/software/swave/Swave-main` | Swave软件路径 |
| `-o, --output-dir` | `./swave_output` | 输出目录路径 |

#### 输入文件选项 | Input File Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--decomposed-vcf` | `None` | Decomposed VCF文件（cactus/pggb必需） |
| `--output-mode` | `auto` | 输出模式 (auto/population/single) |
| `--spec-samples` | `None` | 指定样本列表 |

#### SV检测参数 | SV Detection Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--min-sv-size` | `50` | 最小SV大小（bp） |
| `--max-sv-size` | `1000000` | 最大SV大小（bp） |
| `--max-sv-comps` | `5` | 最大SV组件数 |

#### 处理选项 | Processing Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--dup-to-ins` | `False` | 将duplication报告为insertion |
| `--remove-small` | `False` | 移除小于min_sv_size的节点 |
| `--force-reverse` | `False` | 强制调用反向映射snarls |

#### 性能参数 | Performance Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `1` | 线程数 |

#### 外部工具路径 | External Tool Paths

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--minigraph-path` | `minigraph` | minigraph工具路径 |
| `--gfatools-path` | `gfatools` | gfatools工具路径 |

#### 高级选项 | Advanced Options

| 参数 | 描述 |
|------|------|
| `--spec-snarl` | 只调用特定snarl (例如: '>s1>s2') |
| `--spec-path` | 只调用特定path (例如: '>s3>s4') |

### 其他子命令参数

#### convert-seq 参数

| 参数 | 描述 |
|------|------|
| `--vcf-path` | VCF文件路径（必需） |
| `--gfa-path` | GFA文件路径（必需） |
| `--ref-path` | 参考基因组FASTA文件（必需） |
| `--swave-path` | Swave软件路径 |
| `--output-path` | 输出路径 |
| `--force-pangenie` | 强制输出满足pangenie要求的序列 |

#### convert-plines 参数

| 参数 | 描述 |
|------|------|
| `--gfa-path` | GFA文件路径（必需） |
| `--vcf-path` | VCF文件路径（必需） |
| `--ref-vcf-path` | 参考VCF文件路径 |
| `--swave-path` | Swave软件路径 |
| `--output-path` | 输出路径 |
| `--force-vg` | 强制输出满足vg要求的序列 |

#### extract-csv 参数

| 参数 | 描述 |
|------|------|
| `--vcf-path` | VCF文件路径（必需） |
| `--swave-path` | Swave软件路径 |
| `--spec-csv` | 特定CSV类型（INV/DUP/All） |
| `--output-path` | 输出路径 |

#### extract-sample 参数

| 参数 | 描述 |
|------|------|
| `--vcf-path` | VCF文件路径（必需） |
| `--spec-samples` | 指定样本（必需，可多个） |
| `--swave-path` | Swave软件路径 |
| `--output-path` | 输出路径 |

## 输入文件格式 | Input File Formats

### Assemblies TSV文件 | Assemblies TSV File

样本组装信息文件，每行代表一个样本：

```tsv
NAME    hap1.fasta    hap2.fasta
sample1    /path/to/sample1_hap1.fasta    /path/to/sample1_hap2.fasta
sample2    /path/to/sample2_hap1.fasta    /path/to/sample2_hap2.fasta
```

**文件要求**:
- 第一列为样本名称
- 后续列为该样本的单倍型FASTA文件路径
- 支持多倍体样本（多个单倍型文件）

### GFA文件 | GFA File

泛基因组图文件，可以是：
- **Minigraph格式**: 使用 minigraph 构建
- **Cactus格式**: 使用 minigraph-cactus 构建
- **PGGB格式**: 使用 pggb 构建

### Decomposed VCF文件 | Decomposed VCF File

cactus和pggb图所需的decomposed VCF文件（使用 vg deconstruct 生成）。

## 输出结果 | Output Results

### call 输出文件结构 | call Output File Structure

```
swave_output/
├── swave.log                           # 运行日志
├── swave.sample_level.vcf              # 多等位基因SV
├── swave.sample_level.split.vcf        # 二等位基因SV（拆分后）
└── 99_logs/
    └── swave_pipeline.log              # 完整日志
```

### 输出文件说明 | Output File Description

| 文件 | 描述 |
|------|------|
| `swave.sample_level.vcf` | 多等位基因输出，包含所有检测到的SV |
| `swave.sample_level.split.vcf` | 二等位基因输出，将多等位基因拆分 |
| `swave.log` | Swave工具运行日志 |

### SV质量分级 | SV Quality Classification

Swave提供三级质量标签（QUAL列）：

| 级别 | 描述 | 推荐用途 |
|------|------|----------|
| **HighQual (PASS)** | 高质量SV，预测准确 | 推荐用于下游分析 |
| **MediumQual** | 中等质量，可能存在偏差 | 谨慎使用 |
| **LowQual** | 低质量，可能是假阳性 | 建议过滤 |

## 使用示例 | Usage Examples

### 示例1：完整SV检测流程 | Example 1: Complete SV Detection Pipeline

```bash
# Step 1: 检测SV
biopytools swave call \
    -i human_assemblies.tsv \
    -r t2t.chr20.fa \
    -g pg.minigraph.gfa \
    -s minigraph \
    -o human_sv_detection

# Step 2: 转换序列
biopytools swave convert-seq \
    --vcf-path human_sv_detection/swave.sample_level.vcf \
    --gfa-path pg.minigraph.gfa \
    --ref-path t2t.chr20.fa \
    --output-path human_sv.with_sequences.vcf

# Step 3: 提取特定样本
biopytools swave extract-sample \
    --vcf-path human_sv.with_sequences.vcf \
    --spec-samples sample1 sample2 \
    --output-path subset_samples.vcf
```

### 示例2：Cactus图分析 | Example 2: Cactus Graph Analysis

```bash
# 检测SV
biopytools swave call \
    -i plant_assemblies.tsv \
    -r reference.fa \
    -g pangenome.cactus.gfa \
    -s cactus \
    --decomposed-vcf pangenome.raw.vcf.gz \
    -o plant_sv_detection \
    -t 16

# 提取倒位
biopytools swave extract-csv \
    --vcf-path plant_sv_detection/swave.sample_level.vcf \
    --spec-csv INV \
    --output-path plant_inversions.csv
```

### 示例3：质量控制 | Example 3: Quality Control

```bash
# 严格质量控制
biopytools swave call \
    -i assemblies.tsv \
    -r reference.fa \
    -g pangenome.gfa \
    -s minigraph \
    --min-sv-size 100 \
    --max-sv-size 500000 \
    --max-sv-comps 3 \
    --remove-small \
    -o high_quality_sv

# 只提取HighQual SV
# 过滤LowQual和MediumQual
grep "PASS" high_quality_sv/swave.sample_level.split.vcf > high_quality_sv.filtered.vcf
```

### 示例4：批量处理 | Example 4: Batch Processing

```bash
# 批量检测多个样本
for sample in sample1 sample2 sample3; do
    biopytools swave call \
        -i assemblies.tsv \
        -r reference.fa \
        -g pangenome.gfa \
        -s minigraph \
        --spec-samples $sample \
        -o sv_output_$sample
done
```

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **Swave** (版本 1.0 或更新)
  - 下载地址: https://github.com/songbowang125/Swave
- **Python** (版本 3.7+)
  - PyTorch 1.10.1
  - NumPy
  - pysam
  - opencv-python
  - matplotlib

### 外部工具 | External Tools

- **minigraph** (用于minigraph图)
- **gfatools** (可选)
- **vg** (用于cactus/pggb图的VCF decompose)

### 安装依赖 | Installing Dependencies

```bash
# 创建conda环境
conda env create -f swave_environment.yml

# 激活环境
conda activate swave_v.1.2

# 安装外部工具
conda install -c bioconda minigraph gfatools vg
```

## 注意事项 | Important Notes

1. **GFA来源匹配**: 确保gfa_source参数与实际的GFA文件构建工具一致
2. **Decomposed VCF**: cactus和pggb图必须提供decomposed VCF文件
3. **内存需求**: 大型泛基因组图可能需要大量内存（建议32GB+）
4. **线程配置**: 推荐根据可用CPU核心数设置适当的线程数
5. **质量过滤**: 建议过滤LowQual级别的SV用于下游分析
6. **子命令选择**: 根据分析需求选择合适的子命令

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "GFA文件不存在" 错误**
```bash
# 检查GFA文件路径
ls -lh /path/to/pangenome.gfa

# 确保使用绝对路径或正确的相对路径
biopytools swave call -i assemblies.tsv -r ref.fa -g /full/path/to/graph.gfa -s minigraph
```

**Q: "cactus图需要decomposed_vcf" 错误**
```bash
# 为cactus/pggb图提供decomposed VCF
biopytools swave call ... -s cactus --decomposed-vcf decomposed.vcf.gz
```

**Q: 内存不足错误**
```bash
# 减少线程数
biopytools swave call ... -t 4

# 或使用计算节点提交任务
```

**Q: 检测到大量LowQual SV**
```bash
# 调整参数提高质量
biopytools swave call ... \
    --min-sv-size 100 \
    --max-sv-comps 3 \
    --remove-small
```

**Q: 子命令找不到**
```bash
# 查看所有可用子命令
biopytools swave -h

# 查看特定子命令帮助
biopytools swave call -h
biopytools swave convert-seq -h
```

## 结果解读指南 | Result Interpretation Guide

### VCF格式说明 | VCF Format Description

Swave输出的VCF文件包含以下关键信息：

| 列 | 描述 |
|----|------|
| CHROM | 染色体/序列名称 |
| POS | SV起始位置 |
| ID | SV标识符 |
| REF | 参考序列 |
| ALT | 变异序列 |
| QUAL | 质量分级（PASS/MediumQual/LowQual） |
| FILTER | 过滤状态 |
| INFO | SV信息（类型、大小等） |
| FORMAT | 格式字段 |
| 样本列 | 各样本的基因型 |

### SV类型说明 | SV Type Description

| SV类型 | 描述 | 符号 |
|--------|------|------|
| **DEL** | 缺失 | Deletion |
| **INS** | 插入 | Insertion |
| **DUP** | 重复 | Duplication |
| **INV** | 倒位 | Inversion |
| **BND** | 复杂边界 | Breakend |
| **CPX** | 复杂SV | Complex |

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

**注意**: Swave软件本身为学术非商业使用许可，详见Swave官方说明。

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用Swave相关文献：

```
Swave: Structural variant detection from pangenome graphs
链接: https://github.com/songbowang125/Swave
```
