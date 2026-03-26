# HiFi+Hi-C 基因组组装模块

**基因组组装与去冗余工具 | Genome Assembly and Deduplication Tool**

## 功能概述 | Overview

hifi_hic模块是一个基于HiFi长读长测序数据的基因组组装工具，集成了hifiasm组装和Purge_Dups去冗余功能。该模块支持Hi-C数据辅助的组装质量评估，并提供完整的基因组组装到去冗余的流程。

## 主要特性 | Key Features

- **HiFi组装**: 基于hifiasm的高质量基因组组装
- **去冗余**: 集成Purge_Dups自动去除冗余序列
- **灵活的文件检测**: 自动检测多种Purge_Dups输出格式
- **断点续传**: 支持从任意步骤恢复执行
- **详细日志**: 完整的执行日志和错误追踪
- **双语支持**: 中英文双语注释和日志

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# HiFi组装（默认启用去冗余）
biopytools hifi-hic \
    --hifi hifi_reads.fq.gz \
    -p genome_sample \
    -t 24 \
    -o ./assembly_output

# 禁用去冗余
biopytools hifi-hic \
    --hifi hifi_reads.fq.gz \
    -p genome_sample \
    -t 24 \
    --no-purge-dups \
    -o ./assembly_output
```

### 高级用法 | Advanced Usage

```bash
# 自定义Purge_Dups参数
biopytools hifi-hic \
    --hifi hifi_reads.fq.gz \
    -p genome_sample \
    -t 24 \
    --purge-level 2 \
    --purge-dups-read-type hifi \
    -o ./assembly_output

# 指定Purge_Dups路径
biopytools hifi-hic \
    --hifi hifi_reads.fq.gz \
    -p genome_sample \
    -t 24 \
    --purge-dups-path ~/miniforge3/envs/purge_dups_v.1.2.6 \
    -o ./assembly_output
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `--hifi` | HiFi reads文件 | `--hifi hifi_reads.fq.gz` |
| `-o`, `--output` | 输出目录 | `-o ./output` |

### 组装参数 | Assembly Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-p`, `--prefix` | `genome` | 样本前缀 |
| `-t`, `--threads` | `12` | 线程数 |
| `--genome-size` | `auto` | 基因组大小（例如：100m, 1g） |

### Purge_Dups参数 | Purge_Dups Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--enable-purge-dups` | `True` | 启用去冗余（默认启用） |
| `--purge-level` | `None` | 去冗余级别 (0=不去, 1=轻度, 2/3=深度) |
| `--purge-dups-path` | `~/miniforge3/envs/purge_dups_v.1.2.6` | Purge_Dups软件路径 |
| `--purge-dups-threads` | `同threads` | 去冗余线程数 |
| `--purge-dups-read-type` | `hifi` | reads类型 (pacbio/hifi/illumina) |
| `--hom-cov` | `auto` | 纯合覆盖度 |

### 控制参数 | Control Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--no-purge-dups` | `False` | 禁用去冗余 |
| `--force` | `False` | 强制重新运行 |

## 输出文件 | Output Files

### 目录结构 | Directory Structure

```
assembly_output/
├── 01.raw_output/          # 原始输出
├── 02.fasta/               # 组装结果
│   ├── {prefix}.fa         # 主基因组
│   └── {prefix}.hap.fa     # 单倍型基因组
└── 04.purge_dups/          # 去冗余结果
    └── seqs/               # 去冗余序列
        └── {prefix}_purged.purge.fa
```

### 主要输出文件 | Main Output Files

- **主基因组**: `{prefix}.fa` - 组装的主要基因组序列
- **单倍型**: `{prefix}.hap.fa` - 分离的单倍型序列
- **去冗余基因组**: `{prefix}_purged.purge.fa` - 去除冗余后的高质量基因组

## 工作流程 | Workflow

### 1. HiFi组装 | HiFi Assembly

使用hifiasm进行基因组组装：
- 输入：HiFi reads
- 输出：组装结果（primary.fa和hap.fa）
- 特点：高质量、连续性好

### 2. Purge_Dups去冗余 | Purge_Dups Deduplication

自动去除冗余和重复序列：
- 输入：组装基因组
- 输出：去冗余基因组
- 支持的输出格式：
  - `{prefix}_purged.purge.fa`
  - `{prefix}_purged.hap.fa`

**文件名自动检测**：
- 代码会自动检测多种可能的文件名格式
- 支持的输出目录：`sequences`或`seqs`
- 确保兼容不同版本的Purge_Dups输出

## 示例 | Examples

### 示例1：标准组装流程

```bash
biopytools hifi-hic \
    --hifi /data/hifi_reads.fq.gz \
    -p EcA \
    -t 24 \
    -o ./EcA_assembly
```

**输出**：
- `./EcA_assembly/02.fasta/EcA.fa`
- `./EcA_assembly/04.purge_dups/seqs/EcA_purged.purge.fa`

### 示例2：禁用去冗余

```bash
biopytools hifi-hic \
    --hifi /data/hifi_reads.fq.gz \
    -p EcA \
    -t 24 \
    --no-purge-dups \
    -o ./EcA_assembly
```

**输出**：
- `./EcA_assembly/02.fasta/EcA.fa`
- `./EcA_assembly/02.fasta/EcA.hap.fa`

## 注意事项 | Notes

1. **内存要求**: hifiasm需要较大内存，建议至少100GB RAM
2. **磁盘空间**: 确保有足够的磁盘空间存储中间文件
3. **文件格式**: HiFi reads支持FASTQ和FASTQ.GZ格式
4. **去冗余**: 默认启用去冗余，可使用`--no-purge-dups`禁用
5. **断点续传**: 流程支持断点续传，可随时恢复执行

## 故障排除 | Troubleshooting

### 问题1：找不到去冗余输出文件

**错误信息**：
```
去冗余输出文件未找到|Purged output file not found
```

**解决方案**：
- 检查`04.purge_dups/seqs/`或`04.purge_dups/sequences/`目录
- 代码会自动检测多种文件名格式
- 如果问题仍然存在，请检查Purge_Dups是否成功运行

### 问题2：hifiasm运行失败

**可能原因**：
- 内存不足
- HiFi reads格式错误
- hifiasm未正确安装

**解决方案**：
- 增加内存或减少线程数
- 验证输入文件格式
- 检查hifiasm安装：`which hifiasm`

## 版本历史 | Version History

- **v1.1.0** (2026-03-26)
  - 修复Purge_Dups输出文件名和目录名检测
  - 支持多种文件名格式自动检测
  - 兼容sequences和seqs两种目录名

- **v1.0.0** (2025-xx-xx)
  - 初始版本
  - 集成hifiasm和Purge_Dups
  - 支持断点续传
