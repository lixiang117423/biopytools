# ANNOVAR结果处理功能 | ANNOVAR Results Processing Features

## 概述 | Overview

本文档描述了ANNOVAR模块中新增的结果处理功能，这些功能整合了原有的两个脚本：
- `41.处理Annovar外显子注释结果.py`
- `42.处理Annovar所有注释结果.py`

This document describes the newly added results processing features in the ANNOVAR module, which integrates two original scripts:
- `41.处理Annovar外显子注释结果.py`
- `42.处理Annovar所有注释结果.py`

## 新增的类 | New Classes

### 1. ExonicVariantProcessor
处理外显子变异注释结果的专用处理器 | Specialized processor for exonic variant annotation results.

### 2. AllVariantProcessor
处理所有变异注释结果的专用处理器 | Specialized processor for all variant annotation results.

### 3. ANNOVARResultsProcessor
主要的结果处理器，整合了上述两个功能 | Main results processor that integrates the above two functionalities.

## 使用方法 | Usage Methods

### 方法1：自动处理（推荐）| Method 1: Automatic Processing (Recommended)

当运行完整的注释流程时，结果处理会自动执行：

```python
from biopytools.annovar import ANNOVARAnnotator

# 创建注释器 | Create annotator
annotator = ANNOVARAnnotator(
    gff3_file="annotation.gff3",
    genome_file="genome.fa",
    vcf_file="variants.vcf",
    build_ver="OV"
)

# 运行完整流程，结果会自动处理 | Run full pipeline, results will be automatically processed
annotator.run_full_pipeline()
```

### 方法2：手动处理特定结果 | Method 2: Manual Processing of Specific Results

```python
from biopytools.annovar import ANNOVARAnnotator

annotator = ANNOVARAnnotator(
    gff3_file="annotation.gff3",
    genome_file="genome.fa",
    vcf_file="variants.vcf",
    build_ver="OV"
)

# 仅处理外显子注释结果 | Process exonic annotation results only
exonic_file = annotator.process_exonic_results_only()

# 处理所有变异注释结果 | Process all variant annotation results
all_file = annotator.process_all_results_only()

# 处理所有可用的结果文件 | Process all available result files
processed_files = annotator.process_annotation_results()
```

### 方法3：使用独立的结果处理器 | Method 3: Using Independent Results Processors

```python
from biopytools.annovar import ANNOVARResultsProcessor
import logging

# 设置日志 | Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# 创建结果处理器 | Create results processor
processor = ANNOVARResultsProcessor(logger, "./output")

# 处理外显子注释结果 | Process exonic annotation results
exonic_output = processor.process_exonic_results(
    "path/to/exonic_variant_function",
    "output_prefix"
)

# 处理所有变异注释结果 | Process all variant annotation results
all_output = processor.process_all_results(
    "path/to/variant_function",
    "output_prefix",
    apply_filters=True  # 可选：应用过滤器 | Optional: apply filters
)

# 自动处理所有可用的结果 | Automatically process all available results
processed_files = processor.process_available_results(
    "vcf_basename",
    apply_filters=False
)
```

## 输出文件格式 | Output File Formats

### 外显子变异结果 | Exonic Variant Results
处理后的外显子变异文件包含以下列：
- Line_ID
- 染色体 | Chromosome
- 变异起始 | Variant Start
- 变异终止 | Variant End
- 突变类型 | Mutation Type
- 基因 | Gene
- 转录本 | Transcript
- 变异结果 | Variant Effect
- DNA位置起 | DNA Position From
- DNA位置止 | DNA Position To
- DNA参考 | DNA Reference
- DNA变异 | DNA Variant
- 蛋白位置 | Protein Position
- 蛋白参考 | Protein Reference
- 蛋白变异 | Protein Variant

### 所有变异结果 | All Variant Results
处理后的所有变异文件包含以下列：
- 染色体 | Chromosome
- 起始位置 | Start Position
- 终止位置 | End Position
- 区域类型 | Region Type
- 基因 | Gene
- 距离 | Distance
- 突变类型 | Mutation Type
- 参考序列 | Reference Sequence
- 变异序列 | Variant Sequence
- 变异长度 | Variant Length
- INDEL大小 | INDEL Size
- 频率 | Frequency
- 质量分数 | Quality Score
- 测序深度 | Sequencing Depth

## 过滤选项 | Filtering Options

处理所有变异结果时，可以应用过滤器：

```python
filters = {
    'region_type': ['exonic', 'splicing'],  # 按区域类型过滤 | Filter by region type
    'mutation_type': ['SNP'],               # 按突变类型过滤 | Filter by mutation type
    'exclude_intergenic': True,              # 排除基因间区 | Exclude intergenic regions
    'min_freq': 0.05                         # 最小频率阈值 | Minimum frequency threshold
}

processed_file = processor.process_all_results(
    input_file="variant_function",
    output_prefix="filtered_results",
    apply_filters=True,
    filters=filters
)
```

## 统计信息 | Statistics

处理结果时会自动输出统计信息，包括：
- 区域类型分布 | Region type distribution
- 突变类型分布 | Mutation type distribution
- 基因相关变异数量 | Number of gene-related variants

## 向后兼容性 | Backward Compatibility

所有新的结果处理功能都不会影响现有的ANNOVAR注释流程。原有的方法和参数保持不变，新功能作为可选的增强功能添加。

All new result processing features do not affect the existing ANNOVAR annotation workflow. Original methods and parameters remain unchanged, with new features added as optional enhancements.