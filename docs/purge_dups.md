# Purge_Dups 基因组去冗余分析模块

**基于测序深度的基因组单倍型和重叠序列去除工具 | Genome Haplotigs and Overlaps Removal Tool Based on Read Depth**

## 功能概述 | Overview

Purge_Dups 去冗余分析模块是一个强大的基因组质量控制工具，基于测序深度原理，能够有效识别并去除基因组组装中的单倍型序列（haplotigs）和重叠区域。该工具适用于三代测序（PacBio/HiFi）组装的基因组质量优化，通过去除冗余序列提高基因组的完整性和准确性。

## 主要特性 | Key Features

- **智能深度分析**: 基于测序深度自动计算覆盖度阈值
- **完整去冗余流程**: 比对、深度统计、阈值计算、序列识别、序列提取五步骤自动化
- **灵活步骤控制**: 支持单独运行任意步骤或完整流程
- **多种数据支持**: 支持PacBio CLR、HiFi和Illumina数据
- **详细日志记录**: 完整的处理过程日志和错误追踪
- **高度可配置**: 丰富的参数选项满足不同分析需求

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 使用PacBio数据对基因组进行去冗余分析
biopytools purge-dups \
    -i assembly.fa \
    -r pacbio_reads.fq \
    -o purged_output
```

### 高级用法 | Advanced Usage

```bash
# 使用HiFi数据并自定义参数（默认类型，可省略--read-type）
biopytools purge-dups \
    -i genome.fa \
    -r hifi_reads.fq \
    -t 32 \
    --min-fraction 0.85 \
    -o purged_genome

# 使用PacBio CLR数据（需指定--read-type pacbio）
biopytools purge-dups \
    -i genome.fa \
    -r pacbio_reads.fq \
    --read-type pacbio \
    -t 32 \
    -o purged_genome
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 基因组组装文件(FASTA格式) | `-i assembly.fa` |
| `-r, --reads` | 测序reads文件 | `-r pacbio_reads.fq` |

### 通用配置 | General Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./purge_dups_output` | 输出目录路径 |
| `-t, --threads` | `12` | 线程数 |

### 数据类型配置 | Data Type Configuration

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--read-type` | `hifi` | 测序数据类型 (pacbio/hifi/illumina) |

### 去冗余参数 | Deduplication Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--min-fraction` | `0.8` | 最小比例阈值 (0-1) |
| `--two-round-chaining` | `TRUE` | 启用两轮链式匹配 |
| `--ends-only` | `TRUE` | 只去除contig末端的冗余 (推荐) |
| `--no-ends-only` | - | 也去除contig中间的冗余 (谨慎使用) |
| `--min-primary-length` | `10000` | 最小主contig长度 (bp) |

### 步骤控制 | Step Control

| 参数 | 描述 |
|------|------|
| `-s, --step` | 只运行指定步骤 (1/2/3/4/5) |

### 步骤说明 | Step Descriptions

| 步骤 | 名称 | 描述 |
|------|------|------|
| **1** | 计算深度 | 将reads比对到基因组并计算测序深度 |
| **2** | 计算阈值 | 根据深度分布计算覆盖度阈值 |
| **3** | 分割比对 | 分割基因组并进行自比对 |
| **4** | 去冗余 | 识别冗余序列并生成BED文件 |
| **5** | 获取序列 | 提取纯化后的主序列和单倍型序列 |

### 其他参数 | Other Parameters

| 参数 | 描述 |
|------|------|
| `--split-by-n` | split_fa按N分割scaffold |
| `--manual-cutoffs` | 手动指定阈值文件 |

## 分析流程 | Analysis Pipeline

### 完整流程图 | Complete Pipeline Flow

```
输入: 基因组 + Reads
   ↓
步骤1: Reads比对 → PAF文件 → 深度统计
   ↓
步骤2: 深度分布分析 → 覆盖度阈值
   ↓
步骤3: 基因组分割 → 自比对 → PAF文件
   ↓
步骤4: 冗余识别 → dups.bed文件
   ↓
步骤5: 序列提取 → 主序列 + 单倍型序列
   ↓
输出: 纯化后的基因组
```

## 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
purge_dups_output/
├── coverage/              # 深度统计结果
│   ├── assembly.base.cov  # 碱基级深度文件
│   ├── assembly.stat      # 深度统计文件
│   └── cutoffs            # 覆盖度阈值文件
├── split_aln/             # 分割和比对结果
│   ├── assembly.split     # 分割后的基因组
│   └── assembly.split.self.paf.gz  # 自比对结果
├── purge_dups/            # 去冗余结果
│   ├── dups.bed           # 冗余序列注释文件
│   └── purge_dups.log     # 去冗余日志
├── seqs/                  # 最终序列
│   ├── assembly_purged.purge.fa    # 纯化后的主序列
│   └── assembly_purged.haplotigs.fa  # 单倍型序列
└── purge_dups.log         # 完整日志文件
```

### 关键输出文件说明 | Key Output Files Description

- **`*.purge.fa`**: 去除冗余后的主contig序列，用于后续分析
- **`*.haplotigs.fa`**: 被识别为单倍型的序列，可用于辅助分析
- **`dups.bed`**: 冗余区域的BED格式注释文件
- **`cutoffs`**: 自动计算的覆盖度阈值（低、中、高）

## 使用示例 | Usage Examples

### 示例1: 标准PacBio数据去冗余 | Example 1: Standard PacBio Deduplication

```bash
biopytools purge-dups \
    -i primary_asm.fa \
    -r pacbio_reads.subreads.fq \
    -o pacbio_purged
```

### 示例2: HiFi数据高质量去冗余 | Example 2: HiFi High-Quality Deduplication

```bash
biopytools purge-dups \
    -i hifi_asm.fa \
    -r hifi_reads.fq \
    --min-fraction 0.9 \
    -t 48 \
    -o hifi_purged
```

### 示例3: 分步执行流程 | Example 3: Step-by-Step Execution

```bash
# 第一步: 计算测序深度
biopytools purge-dups -i asm.fa -r reads.fq -s 1 -o step1

# 第二步: 计算覆盖度阈值
biopytools purge-dups -i asm.fa -r reads.fq -s 2 -o step2

# 第三步: 分割基因组并自比对
biopytools purge-dups -i asm.fa -r reads.fq -s 3 -o step3

# 第四步: 识别冗余序列
biopytools purge-dups -i asm.fa -r reads.fq -s 4 -o step4

# 第五步: 获取纯化后的序列
biopytools purge-dups -i asm.fa -r reads.fq -s 5 -o step5
```

### 示例4: 使用手动指定阈值 | Example 4: Using Manual Cutoffs

```bash
# 首先查看深度直方图确定阈值
python3 hist_plot.py -c cutoffs assembly.stat coverage.png

# 使用手动指定的阈值运行
biopytools purge-dups \
    -i asm.fa \
    -r reads.fq \
    --manual-cutoffs my_cutoffs \
    -o manual_purged
```

## 结果解读 | Result Interpretation

### 覆盖度阈值说明 | Coverage Cutoffs Explanation

自动计算的阈值包含三个值：

| 阈值 | 用途 | 说明 |
|------|------|------|
| **低** | JUNK过滤 | 深度低于此值的contig被标记为低质量 |
| **中** | 单倍型检测 | 深度介于低和中值之间的contig检测单倍型 |
| **高** | 重复序列 | 深度高于此值的contig被标记为重复序列 |

### 输出序列分类 | Output Sequence Classification

| 文件 | 内容 | 用途 |
|------|------|------|
| `*.purge.fa` | 纯化后的主序列 | 主要基因组，用于下游分析 |
| `*.haplotigs.fa` | 单倍型序列 | 备用单倍型，可用于特定分析 |

## 注意事项 | Important Notes

1. **推荐参数**:
   - 初次使用建议保留 `--ends-only` 参数，只去除末端冗余
   - 如果怀疑过度去冗余，可以适当提高 `--min-fraction` 值

2. **Illumina数据**:
   - Illumina reads需要先用bwa比对生成BAM文件
   - 建议使用高深度数据（>50x）以获得准确结果

3. **基因组质量**:
   - 输入基因组应为主要组装结果（primary assembly）
   - 避免使用已经经过过滤的基因组

4. **资源需求**:
   - 大基因组（>1Gb）建议使用32+线程和100GB+内存
   - minimap2比对步骤占用内存最多

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: 去冗余后基因组过小**
```bash
# 检查cutoffs阈值是否合理
cat coverage/cutoffs

# 使用手动阈值重新运行
biopytools purge-dups ... --manual-cutoffs manual_cutoffs
```

**Q: 没有去除任何冗余**
```bash
# 降低min_fraction值
biopytools purge-dups ... --min-fraction 0.7
```

**Q: 内存不足**
```bash
# 减少线程数
biopytools purge-dups ... -t 12

# 或分步运行
biopytools purge-dups ... -s 1
# ...继续其他步骤
```

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **Purge_Dups** (版本 1.2.6)
- **Minimap2** (版本 2.24+)
- **Python** (版本 3.7+)
- **Python包**:
  - `click` - 命令行界面

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐24核以上）
- **RAM**: 最少32GB（大基因组推荐100GB以上）
- **存储**: 预留基因组文件大小10倍的磁盘空间

## 相关资源 | Related Resources

- [Purge_Dups GitHub](https://github.com/dfguan/purge_dups)
- [Purge_Dups 原始论文](https://doi.org/10.1093/bioinformatics/btaa353)
- [基因组组装最佳实践](https://github.com/malonge/BlobToolKit)

## 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

## 引用信息 | Citation

如果在学术研究中使用此工具，请引用Purge_Dups相关文献：

```
Guan, D., et al. (2020).
purge_dups: a heterozygous and duplicate removal tool for diploid genome assembly.
Bioinformatics, 36(12), 3756-3760.
```
