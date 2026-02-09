# BAM覆盖度统计工具 | BAM Coverage Statistics Tool

## 版本 | Version: 1.0.0
## 作者 | Author: Xiang LI
## 日期 | Date: 2025-12-31

---

## 功能概述 | Function Overview

BAM覆盖度统计工具用于计算BAM文件中指定染色体区域的每个碱基位置的reads覆盖度。支持单个BAM文件或批量处理目录中的多个BAM文件，并可自动合并多个样本的覆盖度数据。

The BAM Coverage Statistics Tool calculates per-base read coverage for specified chromosomal regions in BAM files. It supports processing single BAM files or batch processing multiple BAM files in a directory, with automatic merging of coverage data from multiple samples.

---

## 主要特性 | Main Features

1. **灵活的输入方式 | Flexible Input**
   - 支持单个BAM文件 | Support single BAM file
   - 支持BAM文件目录 | Support BAM files directory
   - 自动索引管理 | Automatic index management

2. **精确的区域控制 | Precise Region Control**
   - 指定染色体和起止位置 | Specify chromosome and start/end positions
   - 支持分析到染色体末端 | Support analysis to chromosome end
   - 1-based坐标系统 | 1-based coordinate system

3. **质量控制选项 | Quality Control Options**
   - MAPQ过滤 | MAPQ filtering
   - 碱基质量过滤 | Base quality filtering
   - 可自定义质量阈值 | Customizable quality thresholds

4. **完整的输出 | Complete Output**
   - 合并的覆盖度矩阵 | Merged coverage matrix
   - 详细的统计摘要 | Detailed summary statistics
   - 易于导入下游分析 | Easy to import for downstream analysis

---

## 安装要求 | Requirements

### 软件依赖 | Software Dependencies

- Python 3.7+
- samtools (必需 | Required)

### Python包依赖 | Python Package Dependencies

- 无特殊依赖 | No special dependencies
- 使用标准库 | Uses standard library only

---

## 使用方法 | Usage

### 命令行接口 | Command Line Interface

#### 基本用法 | Basic Usage

```bash
# 分析单个BAM文件
biopytools bam-cov -i sample.bam -c chr1 -s 1000000 -e 2000000

# 分析目录中所有BAM文件
biopytools bam-cov -i ./bam_files -c chr1 -s 1000000

```

#### 完整参数示例 | Full Parameter Example

```bash
biopytools bam-cov \\
    -i ./bam_files \\
    -c Chr12 \\
    -s 125000000 \\
    -o coverage_results \\
    -p my_analysis \\
    --min-mapq 20 \\
    --min-baseq 30 \\
    --no-summary \\
    -vv
```

### Python API | Python API

```python
from biopytools.bam_coverage_stats import BAMCoverageAnalyzer

# 创建分析器
analyzer = BAMCoverageAnalyzer(
    input="./bam_files",
    chromosome="chr1",
    start=1000000,
    end=2000000,
    output_prefix="coverage_stats",
    min_mapq=20,
    min_baseq=30
)

# 运行分析
analyzer.run_analysis()
```

---

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 说明 | 示例 |
|------|------|------|
| `-i, --input` | 输入路径（BAM文件或包含BAM的目录）| sample.bam 或 ./bam_files |
| `-c, --chromosome` | 染色体名称 | chr1, Chr12, scaffold_123 |
| `-s, --start` | 起始位置（1-based）| 1000000 |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-e, --end` | 染色体末端 | 终止位置 |
| `--min-mapq` | 0 | 最小mapping质量 (0-255) |
| `--min-baseq` | 0 | 最小碱基质量 (0-93) |
| `-o, --output-dir` | ./bam_coverage_stats_output | 输出目录 |
| `-p, --output-prefix` | coverage | 输出文件前缀 |
| `--no-merge` | False | 不合并样本输出 |
| `--no-summary` | False | 不生成统计摘要 |
| `-v, --verbose` | 0 | 增加输出详细程度 |
| `--quiet` | False | 静默模式 |
| `--log-file` | None | 日志文件路径 |
| `--log-level` | INFO | 日志级别 |
| `--dry-run` | False | 试运行模式 |
| `-t, --threads` | 64 | 线程数 |

---

## 输出文件 | Output Files

### 1. 合并覆盖度文件 | Merged Coverage File

**文件名 | Filename**: `{prefix}_merged.txt`

**格式 | Format**:
```
Chr    Pos        Sample1    Sample2    Sample3
chr1   1000000    25          30          28
chr1   1000001    27          32          29
chr1   1000002    26          31          27
```

**列说明 | Column Description**:
- 第1列: 染色体名称 | Column 1: Chromosome name
- 第2列: 碱基位置 | Column 2: Base position
- 第3列及以后: 各样本的覆盖度 | Column 3+: Coverage for each sample

### 2. 统计摘要文件 | Summary Statistics File

**文件名 | Filename**: `{prefix}_summary.txt`

**格式 | Format**:
```
Sample    Total_Positions    Mean_Coverage    Median_Coverage    ...
Sample1   1000000            25.34            24.00              ...
Sample2   1000000            30.12            29.00              ...
```

**统计指标 | Statistics**:
- Total_Positions: 总位置数
- Mean_Coverage: 平均覆盖度
- Median_Coverage: 中位数覆盖度
- Max_Coverage: 最大覆盖度
- Min_Coverage: 最小覆盖度
- Positions_0X: 覆盖度为0的位点数
- Positions_>0X: 覆盖度>0的位点数
- Positions_>10X: 覆盖度>10的位点数
- Positions_>30X: 覆盖度>30的位点数
- Coverage_%_>0X: 覆盖度>0的位点百分比
- Coverage_%_>10X: 覆盖度>10的位点百分比

### 3. 临时文件 | Temporary Files

**位置 | Location**: `{output_dir}/temp/`

包含每个样本的单独覆盖度文件 | Contains individual coverage file for each sample

---

## 使用场景 | Use Cases

### 1. BSA分析 | BSA Analysis

比较两个池的覆盖度差异，定位目标区域：

```bash
# 分析突变池
biopytools bam-cov -c chr1 -s 1000000 -e 2000000 \\
    -d mutant_pool_bam -o mutant_coverage -p mutant

# 分析野生池
biopytools bam-cov -c chr1 -s 1000000 -e 2000000 \\
    -d wild_pool_bam -o wild_coverage -p wild
```

### 2. 重测序数据质控 | Resequencing Data QC

评估测序覆盖度：

```bash
biopytools bam-cov -i sample.bam -c chr1 -s 1 -e 50000000 \\
    --min-mapq 20 --min-baseq 30
```

### 3. 目标区域覆盖度分析 | Target Region Coverage

分析特定基因或区域的覆盖度：

```bash
biopytools bam-cov -c chr5 -s 15000000 -e 15050000 \\
    -d ./samples_bam -p gene_xyz_coverage
```

---

## 质量过滤建议 | Quality Filtering Recommendations

### MAPQ (Mapping Quality)

| 值 | 用途 | 说明 |
|----|------|------|
| 0 | 保留所有比对 | Keep all alignments |
| 10 | 宽松过滤 | Loose filtering |
| 20 | 推荐用于变异检测 | Recommended for variant calling |
| 30 | 高质量分析 | High-quality analysis |

### Base Quality

| 值 | 用途 | 说明 |
|----|------|------|
| 0 | 保留所有碱基 | Keep all bases |
| 20 | 标准质量过滤 | Standard quality filtering |
| 30 | 严格质量过滤 | Strict quality filtering |

---

## 性能优化 | Performance Optimization

### 1. BAM索引 | BAM Indexing

- 确保BAM文件已索引 | Ensure BAM files are indexed
- 自动创建缺失的索引 | Automatically create missing indexes
- 索引文件格式：.bai或.bami | Index file format: .bai or .bami

### 2. 区域大小 | Region Size

- 小区域（<1Mb）: 快速处理 | Small regions (<1Mb): Fast processing
- 中等区域（1-10Mb）: 中等速度 | Medium regions (1-10Mb): Moderate speed
- 大区域（>10Mb）: 需要更多内存 | Large regions (>10Mb): More memory needed

### 3. 批量处理 | Batch Processing

- 使用`-i`参数指定目录批量处理 | Use `-i` parameter to specify directory for batch processing
- 自动合并多个样本 | Automatically merge multiple samples
- 生成统一的输出矩阵 | Generate unified output matrix

---

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**1. BAM文件未索引 | BAM file not indexed**
```
解决 | Solution: 程序会自动创建索引
The program will automatically create the index
```

**2. 染色体名称不匹配 | Chromosome name mismatch**
```
解决 | Solution: 检查BAM头部的染色体名称
Solution: Check chromosome name in BAM header
samtools view -H file.bam | grep @SQ
```

**3. 内存不足 | Out of memory**
```
解决 | Solution: 减小分析区域或提高过滤阈值
Solution: Reduce analysis region or increase filtering threshold
```

**4. samtools未找到 | samtools not found**
```
解决 | Solution: 确保samtools在PATH中
Solution: Ensure samtools is in PATH
which samtools
```

---

## 示例工作流 | Example Workflow

### BSA分析工作流 | BSA Analysis Workflow

```bash
#!/bin/bash

# 步骤1: 分析突变池
echo "分析突变池 | Analyzing mutant pool"
biopytools bam-cov \\
    -c Chr12 \\
    -s 125000000 \\
    -d ./mutant_bam \\
    -o bs_a_results \\
    -p mutant \\
    --min-mapq 20 \\
    --min-baseq 30

# 步骤2: 分析野生池
echo "分析野生池 | Analyzing wild pool"
biopytools bam-cov \\
    -c Chr12 \\
    -s 125000000 \\
    -d ./wild_bam \\
    -o bs_a_results \\
    -p wild \\
    --min-mapq 20 \\
    --min-baseq 30

# 步骤3: 计算覆盖度比值比（使用R或其他工具）
echo "计算覆盖度比值 | Calculating coverage ratio"
# 使用R脚本计算两个池之间的覆盖度差异
# Use R script to calculate coverage difference between pools
```

---

## 版本历史 | Version History

### v1.0.0 (2025-12-31)
- 初始版本 | Initial release
- 支持单个和批量BAM文件处理 | Support single and batch BAM processing
- 覆盖度统计和合并功能 | Coverage statistics and merging functionality
- 详细统计摘要生成 | Detailed summary statistics generation

---

## 联系方式 | Contact

作者 | Author: Xiang LI
项目 | Project: biopytools
日期 | Date: 2025-12-31
