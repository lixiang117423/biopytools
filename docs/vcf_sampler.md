# VCF抽样工具 | VCF Sampling Tool

## 功能描述 | Function Description

从VCF文件中按比例随机抽取SNP位点。该工具针对每条染色体分别进行抽样，确保每条染色体的抽样比例一致。

Randomly sample SNP sites from VCF files by proportion. This tool performs sampling separately for each chromosome to ensure consistent sampling rates across chromosomes.

## 版本信息 | Version

- **版本 | Version**: 1.0.0
- **作者 | Author**: Xiang LI
- **日期 | Date**: 2025-01-03

## 模块结构 | Module Structure

```
vcf_sampler/
├── __init__.py       # 模块入口 | Module entry
├── config.py         # 配置管理 | Configuration management
├── utils.py          # 工具函数 | Utility functions
├── sampler.py        # 核心抽样逻辑 | Core sampling logic
└── main.py           # 主程序 | Main program
```

## 使用方法 | Usage

### 命令行使用 | Command Line Usage

```bash
# 基本用法 - 抽取25%的SNP (使用默认随机种子1288) | Basic usage - Sample 25% of SNPs (using default random seed 1288)
biopytools vcf-sampler -i input.vcf.gz -o output.vcf.gz

# 抽取5%的SNP | Sample 5% of SNPs
biopytools vcf-sampler -i input.vcf.gz -o output.vcf.gz -r 0.05

# 使用自定义随机种子确保可重复性 | Use custom random seed for reproducibility
biopytools vcf-sampler -i input.vcf.gz -o output.vcf.gz -r 0.1 -s 12345

# 抽取50%的SNP并保存日志 | Sample 50% of SNPs and save log
biopytools vcf-sampler -i input.vcf.gz -o output.vcf.gz -r 0.5 --log-file sample.log -v
```

### Python API使用 | Python API Usage

```python
from biopytools.vcf_sampler import VCFSampler, VCFSamplerConfig

# 方式1: 直接创建抽样器 | Method 1: Create sampler directly
sampler = VCFSampler(
    input_vcf="input.vcf.gz",
    output_vcf="output.vcf.gz",
    sample_rate=0.25,  # 抽样比例 | Sampling rate
    random_seed=1288  # 随机种子 (默认1288) | Random seed (default 1288)
)
sampler.run_sampling()

# 方式2: 使用配置对象 | Method 2: Use configuration object
config = VCFSamplerConfig(
    input_vcf="input.vcf.gz",
    output_vcf="output.vcf.gz",
    sample_rate=0.25,
    random_seed=1288
)
sampler = VCFSampler(config=config)
sampler.run_sampling()
```

## 参数说明 | Parameters

| 参数 | Parameter | 说明 | Description | 默认值 | Default |
|------|-----------|------|-------------|--------|---------|
| `-i, --input` | 输入VCF文件 | Input VCF file (必须是.vcf.gz格式 | must be .vcf.gz) | 必需 | Required |
| `-o, --output` | 输出VCF文件 | Output VCF file (必须是.vcf.gz格式 | must be .vcf.gz) | 必需 | Required |
| `-r, --sample-rate` | 抽样比例 | Sampling rate (0.0-1.0) | 0.25 | 0.25 |
| `-s, --random-seed` | 随机种子 | Random seed for reproducibility | 1288 | 1288 |
| `--log-file` | 日志文件 | Log file path | None | None |
| `-v, --verbose` | 详细输出 | Verbose mode (-v: INFO, -vv: DEBUG) | False | False |

## 抽样流程 | Sampling Process

该工具采用三遍扫描策略 | This tool uses a three-pass scanning strategy:

1. **第一遍扫描 | First Pass**: 统计每条染色体的SNP数量 | Count SNPs per chromosome
   - 使用pysam时: 通过索引快速统计 (秒级) | With pysam: Fast counting via index (seconds)
   - 使用Python时: 逐行扫描 (分钟级) | With Python: Line-by-line scan (minutes)
2. **第二遍扫描 | Second Pass**: 随机选择要抽取的SNP索引 | Randomly select SNP indices to sample
3. **第三遍扫描 | Third Pass**: 写入抽样的VCF文件 | Write sampled VCF file
   - 使用pysam时: 直接访问特定染色体记录 | With pysam: Direct chromosome access
   - 使用Python时: 逐行过滤写入 | With Python: Line-by-line filtering

### 性能对比 | Performance Comparison

以100万个SNP的VCF文件为例，抽取25%的SNP：

| 模式 | Mode | 统计时间 | Counting | 写入时间 | Writing | 总时间 | Total |
|------|------|---------|----------|---------|----------|--------|--------|
| **pysam** | ~2秒 | ~5秒 | **~7秒** | ~7 seconds |
| **Python** | ~60秒 | ~120秒 | **~180秒** | ~180 seconds |

**性能提升**: 使用pysam可以提升 **10-25倍** 的速度！| **Speedup**: Using pysam can achieve **10-25x** faster processing!

## 输出示例 | Output Example

```
[2025-01-03 10:30:15] INFO: ============================================================
[2025-01-03 10:30:15] INFO: VCF抽样流程 | VCF Sampling Pipeline
[2025-01-03 10:30:15] INFO: ============================================================
[2025-01-03 10:30:15] INFO: 输入文件 | Input file: /path/to/input.vcf.gz
[2025-01-03 10:30:15] INFO: 输出文件 | Output file: /path/to/output.vcf.gz
[2025-01-03 10:30:15] INFO: 抽样比例 | Sampling rate: 25.0%
[2025-01-03 10:30:15] INFO: 随机种子 | Random seed: 1288
[2025-01-03 10:30:15] INFO: ============================================================
[2025-01-03 10:30:15] INFO: 步骤 1/3 | Step 1/3: 统计SNP数量 | Count SNPs by chromosome
[2025-01-03 10:30:15] INFO: ============================================================
[2025-01-03 10:30:15] INFO: VCF索引不存在，开始创建索引 | VCF index not found, creating index...
[2025-01-03 10:30:15] INFO: 文件 | File: /path/to/input.vcf.gz
[2025-01-03 10:30:20] INFO: 索引创建成功 | Index created successfully
[2025-01-03 10:30:20] INFO: 总SNP数量 | Total SNP count: 1,000,000
[2025-01-03 10:30:20] INFO: 染色体数量 | Chromosome count: 20
[2025-01-03 10:30:20] INFO: 步骤 1 完成 | Step 1 completed
[2025-01-03 10:30:20] INFO: 步骤 2/3 | Step 2/3: 选择SNP索引 | Select SNP indices
[2025-01-03 10:30:20] INFO: ============================================================
[2025-01-03 10:30:20] INFO: 开始随机抽取SNP (抽样比例 | sampling rate: 25.0%)
[2025-01-03 10:30:20] INFO: 总共选择 | Total selected: 250,000 个SNP | SNPs
[2025-01-03 10:30:20] INFO: 步骤 2 完成 | Step 2 completed
[2025-01-03 10:30:20] INFO: 步骤 3/3 | Step 3/3: 写入抽样结果 | Write sampled results
[2025-01-03 10:30:20] INFO: ============================================================
[2025-01-03 10:30:25] INFO: 总SNP数 | Total SNP count: 1,000,000
[2025-01-03 10:30:25] INFO: 写入的SNP数 | Written SNP count: 250,000
[2025-01-03 10:30:25] INFO: 实际抽样比例 | Actual sampling rate: 25.00% (目标 | target: 25.0%)
[2025-01-03 10:30:25] INFO: 输出文件 | Output file: /path/to/output.vcf.gz
[2025-01-03 10:30:25] INFO: 创建VCF索引 | Creating VCF index...
[2025-01-03 10:30:28] INFO: 索引创建成功 | Index created successfully
[2025-01-03 10:30:28] INFO: 已删除自动创建的索引 | Removed auto-created index: /path/to/input.vcf.gz.tbi
[2025-01-03 10:30:28] INFO: 步骤 3 完成 | Step 3 completed
[2025-01-03 10:30:28] INFO: ============================================================
[2025-01-03 10:30:28] INFO: 抽样总结 | Sampling Summary
[2025-01-03 10:30:28] INFO: ============================================================
[2025-01-03 10:30:28] INFO: 总SNP数 | Total SNPs: 1,000,000
[2025-01-03 10:30:28] INFO: 抽取的SNP数 | Sampled SNPs: 250,000
[2025-01-03 10:30:28] INFO: 实际抽样比例 | Actual sampling rate: 25.00% (目标 | target: 25.0%)
[2025-01-03 10:30:28] INFO: 运行时间 | Runtime: 13.23 秒 | seconds
[2025-01-03 10:30:28] INFO: 输出文件 | Output file: /path/to/output.vcf.gz
[2025-01-03 10:30:28] INFO: ============================================================
[2025-01-03 10:30:28] INFO: VCF抽样流程成功完成 | VCF sampling pipeline completed successfully
```

## 注意事项 | Notes

1. **文件格式**: 输入和输出文件必须是gzip压缩的VCF格式 (.vcf.gz) | **File Format**: Input and output files must be gzip-compressed VCF format (.vcf.gz)
2. **自动索引**: 如果输入VCF没有索引，工具会自动使用tabix创建，并在完成后自动删除临时索引 | **Auto Index**: Tool will auto-create index using tabix if missing, and auto-cleanup after completion
3. **依赖要求**: 要获得最佳性能，需要安装pysam和htslib (tabix) | **Dependencies**: For best performance, install pysam and htslib (tabix)
4. **抽样比例**: 必须在0.0到1.0之间，默认为0.25 | **Sample Rate**: Must be between 0.0 and 1.0, default is 0.25
5. **随机种子**: 默认使用1288作为随机种子，可自定义以确保结果可重复 | **Random Seed**: Default seed is 1288, can be customized for reproducible results
6. **内存占用**: 需要存储每条染色体的SNP索引集合 | **Memory Usage**: Needs to store SNP index sets for each chromosome
7. **输出目录**: 输出文件的目录会自动创建 | **Output Directory**: Output directory will be created automatically

## 依赖要求 | Requirements

- Python >= 3.6
- gzip (标准库 | Standard library)
- random (标准库 | Standard library)
- collections (标准库 | Standard library)
- **pysam** (可选 | Optional, **强烈推荐** | **Highly Recommended**)

### 性能优化 | Performance Optimization

该工具支持两种运行模式：

1. **pysam模式 (推荐)**: 使用pysam库，利用VCF索引快速统计和访问，**速度提升10-100倍**
   - 安装pysam: `pip install pysam` 或 `conda install -c bioconda pysam`
   - 安装htslib (包含tabix): `conda install -c bioconda htslib`
   - **自动索引**: 如果输入VCF没有索引，工具会自动使用tabix创建
   - **自动清理**: 工具完成后会自动删除临时创建的索引文件
   - 优点: 快速、高效、自动创建输出索引

2. **Python模式 (兼容模式)**: 纯Python实现，无需额外依赖
   - 优点: 无需安装pysam，兼容性好
   - 缺点: 处理大文件时较慢

如果没有安装pysam或无法创建索引，工具会自动回退到Python模式并给出警告。

## 原始脚本 | Original Script

该模块基于以下脚本重构 | This module is refactored from the following script:
- `/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/05.所有样品/08.10%的SNP的admixture/sample_vcf.py`

## 更新日志 | Changelog

### v1.0.0 (2025-01-03)
- 初始版本 | Initial version
- 实现基础VCF抽样功能 | Implemented basic VCF sampling functionality
- 添加命令行接口 | Added command-line interface
- 添加Python API接口 | Added Python API interface
- 添加详细日志输出 | Added detailed logging output
- 添加随机种子支持 | Added random seed support
