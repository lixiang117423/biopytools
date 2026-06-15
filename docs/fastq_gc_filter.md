# FASTQ文件GC含量和序列长度过滤工具

**专业的FASTQ文件质量控制工具 | Professional FASTQ File Quality Control Tool**

## 功能概述 | Overview

FASTQ文件GC含量和序列长度过滤工具是一个高效的序列质量控制工具,能够根据GC含量和序列长度快速过滤FASTQ文件中的reads。支持普通文本和gzip压缩格式,适用于三代基因组组装前的read质量控制,通过过滤GC含量异常或长度不合适的reads来提高组装质量。

## 主要特性 | Key Features

- **双重过滤机制**: 同时支持GC含量和序列长度过滤
- **灵活参数配置**: 自定义GC含量范围和序列长度范围
- **压缩文件支持**: 自动识别并处理gzip压缩的FASTQ文件(.gz)
- **高效处理**: 流式读取处理,支持大规模FASTQ文件
- **详细统计信息**: 提供完整的过滤统计和通过率信息
- **标准日志输出**: 规范化的日志格式,便于追踪处理进度
- **中英文对照**: 所有输出信息均采用中英文对照格式

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 过滤GC含量大于25%且序列长度大于50bp的reads
biopytools fastq-gc-filter \
    -i input.fastq \
    -o filtered.fastq

# 处理gzip压缩文件
biopytools fastq-gc-filter \
    -i input.fastq.gz \
    -o filtered.fastq.gz \
    --min-gc 30 \
    --min-length 100
```

### 高级用法 | Advanced Usage

```bash
# 精确控制GC含量范围和序列长度范围
biopytools fastq-gc-filter \
    -i reads.fastq.gz \
    -o filtered.fastq.gz \
    --min-gc 30 \
    --max-gc 70 \
    --min-length 100 \
    --max-length 30000

# 过滤低GC含量和短序列
biopytools fastq-gc-filter \
    -i raw_reads.fastq \
    -o clean_reads.fastq \
    --min-gc 40 \
    --min-length 200
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入FASTQ文件路径 | `-i input.fastq` 或 `-i input.fastq.gz` |
| `-o, --output` | 输出FASTQ文件路径 | `-o filtered.fastq` 或 `-o filtered.fastq.gz` |

### GC含量过滤参数 | GC Content Filtering Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--min-gc` | `25.0` | 最小GC含量百分比 (0-100) |
| `--max-gc` | `100.0` | 最大GC含量百分比 (0-100) |

### 序列长度过滤参数 | Sequence Length Filtering Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--min-length` | `50` | 最短序列长度 (bp) |
| `--max-length` | `None` | 最长序列长度 (bp), None表示不限制 |

## 使用场景 | Use Cases

### 三代基因组组装质量控制
```bash
# 过滤HiFi reads,保留高质量序列
biopytools fastq-gc-filter \
    -i hifi_reads.fastq.gz \
    -o filtered_hifi.fastq.gz \
    --min-gc 35 \
    --max-gc 65 \
    --min-length 100
```

### 去除短序列和异常GC含量序列
```bash
# 同时过滤短序列和GC含量异常的序列
biopytools fastq-gc-filter \
    -i raw_data.fastq \
    -o clean_data.fastq \
    --min-gc 25 \
    --max-gc 75 \
    --min-length 200 \
    --max-length 50000
```

## 输出示例 | Output Example

```
2026-02-03 12:00:00.123 - INFO - 开始过滤FASTQ文件|Starting FASTQ filtering
2026-02-03 12:00:00.124 - INFO - 输入文件|Input file: /path/to/input.fastq
2026-02-03 12:00:00.124 - INFO - 输出文件|Output file: /path/to/output.fastq
2026-02-03 12:00:00.124 - INFO - GC含量范围|GC content range: 25.0% - 100.0%
2026-02-03 12:00:00.124 - INFO - 序列长度范围|Sequence length range: 50 - None
2026-02-03 12:00:10.456 - INFO - --------------------------------------------------
2026-02-03 12:00:10.456 - INFO - 处理完成|Processing completed
2026-02-03 12:00:10.456 - INFO - 总reads数|Total reads: 1.50M
2026-02-03 12:00:10.456 - INFO - 通过筛选|Passed filters: 1.45M
2026-02-03 12:00:10.456 - INFO - 过滤掉|Filtered out: 50.00K
2026-02-03 12:00:10.456 - INFO - 通过率|Pass rate: 96.67%
2026-02-03 12:00:10.456 - INFO - 输出文件|Output file: /path/to/output.fastq
```

## 技术细节 | Technical Details

### GC含量计算
- 序列转换为大写后统计G和C碱基数量
- GC含量 = (G碱基数 + C碱基数) / 序列长度 × 100%

### 序列长度过滤
- 最小长度默认值为50bp,可根据需求调整
- 最大长度默认为None(不限制),适用于超长reads(如PacBio HiFi)

### 文件格式支持
- 自动检测.gz扩展名并使用gzip模块处理
- 保留原始文件的压缩状态(输入压缩则输出压缩)

## 注意事项 | Notes

1. **参数范围**: GC含量必须在0-100之间,序列长度必须为正整数
2. **参数顺序**: 确保最小值不大于最大值
3. **内存使用**: 采用流式读取,内存占用较低,适合处理大文件
4. **文件覆盖**: 输出文件如果已存在将被覆盖

## 版本信息 | Version Information

- **版本|Version**: 1.0.0
- **作者|Author**: Xiang LI
- **日期|Date**: 2026-02-03
- **依赖|Dependencies**: Python 3.6+, gzip, logging, pathlib

## 相关工具 | Related Tools

- `biopytools fastp` - FASTQ数据质量控制
- `biopytools fq-stats` - FASTQ文件统计工具
- `biopytools bam2fastq` - BAM to FASTQ转换
