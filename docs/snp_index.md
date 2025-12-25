# SNP Index计算和分析工具 | SNP Index Calculation and Analysis Tool

版本 | Version: 1.0.0
作者 | Author: Xiang LI
日期 | Date: 2025-12-20

## 概述 | Overview

SNP Index工具是一个专门用于**BSA（Bulked Segregant Analysis）**分析的Python包，可以从VCF文件计算SNP index和ΔSNP index，提供完整的统计分析、目标区域检测和可视化功能。

The SNP Index tool is a Python package specifically designed for **BSA (Bulked Segregant Analysis)** that calculates SNP index and ΔSNP index from VCF files, providing complete statistical analysis, target region detection, and visualization functionality.

## 功能特点 | Features

- 🧮 **SNP Index计算**: 从VCF文件计算SNP index和ΔSNP index
- 📊 **统计分析**: 详细的统计信息和分布分析
- 🗺️ **区域检测**: 自动检测潜在的目标区域
- 📈 **可视化**: 多种类型的图表和曼哈顿图
- 🎯 **滑动窗口分析**: 支持滑动窗口折线图（默认启用）
- 📊 **置信区间**: 自动计算和显示置信区间
- 🔧 **灵活过滤**: 支持深度、质量等多种过滤条件
- 📋 **多种输出**: 支持TSV、PNG等多种输出格式
- 🚀 **高性能**: 支持大文件处理和进度显示

## 计算原理 | Calculation Principle

### SNP Index计算公式 | SNP Index Formula

```
SNP index = Alt_depth / (Ref_depth + Alt_depth)
```

其中：
- Alt_depth: 替代等位基因测序深度 | Alternative allele sequencing depth
- Ref_depth: 参考等位基因测序深度 | Reference allele sequencing depth

### ΔSNP Index计算公式 | ΔSNP Index Formula

```
ΔSNP index = SNP_index_pool1 - SNP_index_pool2
```

## 安装和使用 | Installation and Usage

### 作为biopytools模块使用 | Using as biopytools module

```bash
# 完整分析流程
biopytools snp-index -i input.vcf.gz -o results/ -v

# 自定义参数
biopytools snp-index -i input.vcf.gz -o results/ \
    --min-depth 20 --extreme-threshold 0.9

# 只分析已有结果
biopytools snp-index --analyze-only -r results.tsv -o analysis/
```

### 作为Python模块使用 | Using as Python module

```python
from biopytools.snp_index import SNPIndexConfig, SNPIndexProcessor

# 创建配置
config = SNPIndexConfig(
    input_vcf="input.vcf.gz",
    output_dir="results",
    min_depth=10,
    min_quality=20
)

# 运行分析
processor = SNPIndexProcessor(config)
success = processor.run_full_pipeline()
```

## 命令行参数 | Command Line Arguments

### 必需参数 | Required Arguments

| 参数 | 说明 | 示例 |
|------|------|------|
| `-i, --input` | 输入VCF文件路径 | `-i input.vcf.gz` |

### 通用参数 | General Arguments

| 参数 | 默认值 | 说明 | 示例 |
|------|--------|------|------|
| `-o, --output` | `./snp_index_output` | 输出目录 | `-o results/` |
| `-p, --prefix` | `snp_index` | 输出文件前缀 | `-p my_analysis` |
| `-r, --result-file` | None | 已有结果文件 | `-r results.tsv` |

### 过滤参数 | Filtering Parameters

| 参数 | 默认值 | 说明 | 示例 |
|------|--------|------|------|
| `--min-depth` | 10 | 最小测序深度 | `--min-depth 20` |
| `--min-quality` | 20 | 最小质量值 | `--min-quality 30` |
| `--min-mapping-quality` | 20 | 最小mapping质量 | `--min-mapping-quality 30` |
| `--sample-names` | None | 指定样本名称 | `--sample-names pool1 pool2` |

### 分析参数 | Analysis Parameters

| 参数 | 默认值 | 说明 | 示例 |
|------|--------|------|------|
| `--extreme-threshold` | 0.8 | 极端ΔSNP index阈值 | `--extreme-threshold 0.9` |
| `--region-threshold` | 0.5 | 区域检测阈值 | `--region-threshold 0.6` |
| `--min-region-snps` | 5 | 区域最少SNP数量 | `--min-region-snps 10` |
| `--max-region-gap` | 10000 | 区域最大gap(bp) | `--max-region-gap 50000` |

### 模式选择 | Mode Selection

| 参数 | 说明 |
|------|------|
| `--calculate-only` | 只计算SNP index，不进行分析 |
| `--analyze-only` | 只分析已有结果，不计算 |
| `--skip-visualization` | 跳过可视化图表生成 |
| `--disable-sliding-window-plot` | 禁用滑动窗口折线图（默认启用） |

## 输出文件 | Output Files

### 主要结果文件 | Main Result Files

1. **{prefix}_results.tsv**: SNP index计算结果
   - Chromosome, Position, Reference, Alternative
   - Sample1_Ref_Depth, Sample1_Alt_Depth, Sample1_SNP_index
   - Sample2_Ref_Depth, Sample2_Alt_Depth, Sample2_SNP_index
   - Delta_SNP_index

2. **{prefix}_extreme_sites.tsv**: 极端ΔSNP index位点
   - |ΔSNP index| > extreme_threshold的位点

3. **{prefix}_potential_regions.tsv**: 潜在目标区域
   - 连续的极端位点组成的候选区域

4. **{prefix}_sliding_windows.tsv**: 滑动窗口计算结果
   - 每个窗口的统计信息和ΔSNP index
5. **{prefix}_confidence_intervals.txt**: 置信区间分析结果
6. **{prefix}_candidate_regions.tsv**: 候选区域识别结果

### 可视化图表 | Visualization Plots

1. **{prefix}_comprehensive.png**: 综合分析图
   - ΔSNP index分布直方图
   - SNP index分布对比
   - 样本间散点图
   - 染色体分布图

2. **{prefix}_manhattan.png**: 曼哈顿图
   - 全基因组ΔSNP index分布

3. **{prefix}_delta_distribution.png**: ΔSNP index详细分布

4. **{prefix}_correlation.png**: 样本间相关性图

5. **{prefix}_sliding_window.png**: 滑动窗口折线图（默认生成）
   - 全基因组滑动窗口分析
   - 包含置信区间和显著性阈值线
   - 染色体分隔和标记

6. **{prefix}_multi_chrom_sliding.png**: 多染色体分离图 (需要--create-multi-chrom-plot)
   - 每个染色体独立的滑动窗口分析
   - 便于详细查看单个染色体的模式

## 使用示例 | Usage Examples

### 基本使用 | Basic Usage

```bash
# 完整分析流程
biopytools snp-index -i variation.filtered.snp.vcf.gz -o results/ -v
```

### 高级使用 | Advanced Usage

```bash
# 自定义所有参数
biopytools snp-index \
    -i input.vcf.gz \
    -o results/ \
    -p project_x \
    --min-depth 20 \
    --min-quality 30 \
    --extreme-threshold 0.9 \
    --region-threshold 0.6 \
    --min-region-snps 8 \
    --sample-names bulk_A bulk_B \
    -v
```

### 只分析已有结果 | Analyze Existing Results Only

```bash
# 只进行统计分析和可视化
biopytools snp-index \
    --analyze-only \
    -r existing_results.tsv \
    -o analysis_output/ \
    -v
```

### 滑动窗口分析 | Sliding Window Analysis

```bash
# 滑动窗口分析（默认启用）
biopytools snp-index \
    -i input.vcf.gz \
    -o results/ \
    --window-size 500000 \
    --step-size 50000 \
    --confidence-level 0.95 \
    -v

# 创建多染色体分离图
biopytools snp-index \
    -i input.vcf.gz \
    -o results/ \
    --create-multi-chrom-plot \
    --window-size 1000000 \
    --step-size 100000 \
    -v

# 禁用滑动窗口分析
biopytools snp-index \
    -i input.vcf.gz \
    -o results/ \
    --disable-sliding-window-plot \
    -v

# 同时生成多染色体分离图
biopytools snp-index \
    -i input.vcf.gz \
    -o results/ \
    --create-multi-chrom-plot \
    --window-size 1000000 \
    --step-size 100000 \
    --min-window-snps 10 \
    --confidence-level 0.99 \
    -v
```

### 只计算不分析 | Calculate Only

```bash
# 只计算SNP index，跳过分析和可视化
biopytools snp-index \
    -i input.vcf.gz \
    -o results/ \
    --calculate-only
```

## Python API使用 | Python API Usage

### 基本示例 | Basic Example

```python
from biopytools.snp_index import SNPIndexConfig, SNPIndexProcessor

# 创建配置
config = SNPIndexConfig(
    input_vcf="input.vcf.gz",
    output_dir="results",
    prefix="my_analysis",
    min_depth=10,
    min_quality=20
)

# 创建处理器
processor = SNPIndexProcessor(config)

# 运行完整流程
success = processor.run_full_pipeline()
if success:
    print("分析完成！")
```

### 分步骤使用 | Step-by-step Usage

```python
from biopytools.snp_index import (
    SNPIndexConfig, SNPIndexCalculator,
    SNPIndexAnalyzer, SNPIndexVisualizer
)

# 1. 计算SNP index
config = SNPIndexConfig(input_vcf="input.vcf.gz", output_file="results.tsv")
calculator = SNPIndexCalculator(config)
calculator.calculate()

# 2. 分析结果
analyzer = SNPIndexAnalyzer("results.tsv", config)
analyzer.analyze()

# 3. 创建可视化
visualizer = SNPIndexVisualizer(analyzer.data, analyzer.sample_names, config)
visualizer.create_comprehensive_plot("comprehensive.png")
visualizer.create_manhattan_plot("manhattan.png")
```

## 结果解读 | Result Interpretation

### ΔSNP Index解读 | ΔSNP Index Interpretation

- **ΔSNP index ≈ 0**: 两个池在该位点无差异
- **ΔSNP index > 0**: Pool1中Alt等位基因频率更高
- **ΔSNP index < 0**: Pool2中Alt等位基因频率更高

### 极端位点 | Extreme Sites

|ΔSNP index| > 0.8 或 < -0.8 的位点被认为是极端位点，可能与目标性状相关。

### 潜在目标区域 | Potential Target Regions

连续的极端位点（默认5个以上，|ΔSNP index| > 0.5）组成的区域，被认为是潜在的目标区域。

## 性能优化 | Performance Optimization

### 大文件处理 | Large File Processing

- 使用gzip压缩的VCF文件以节省空间
- 调整`--min-depth`参数过滤低质量位点
- 使用`--calculate-only`模式只进行必要计算

### 内存优化 | Memory Optimization

- 工具采用流式处理，支持大文件
- 可视化部分可使用`--skip-visualization`跳过

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

1. **ImportError: matplotlib not found**
   ```bash
   pip install matplotlib
   ```

2. **VCF文件格式错误**
   - 确保VCF文件包含至少2个样本
   - 检查AD字段是否存在

3. **内存不足**
   - 增加过滤条件（如`--min-depth 20`）
   - 使用`--calculate-only`模式

4. **中文字体显示问题**
   - 所有图表现在只使用英文标签，避免字体渲染问题
   - 生成的图中不再显示中文字符

5. **输出目录权限问题**
   ```bash
   mkdir -p output_directory
   chmod 755 output_directory
   ```

## 技术支持 | Technical Support

如有问题或建议，请联系：
For questions or suggestions, please contact:

- 邮箱 | Email: [your-email@example.com]
- 项目主页 | Project: [project-url]

## 更新日志 | Changelog

### v1.1.0 (2025-12-20)
- 🎯 滑动窗口折线图默认启用
- 📊 添加滑动窗口结果保存功能
- 🌐 图表全部使用英文标签，避免字体问题
- 🔧 新增 `--disable-sliding-window-plot` 参数
- 📁 新增输出文件：滑动窗口结果、置信区间、候选区域
- 🐛 修复命名和配置问题

### v1.0.0 (2025-12-20)
- 🎉 初始版本发布
- ✨ 完整的SNP index计算和分析功能
- 📊 多种可视化图表
- 🔧 灵活的参数配置
- 📚 详细的文档和示例