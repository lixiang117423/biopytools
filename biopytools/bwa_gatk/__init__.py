"""
🧬 BWA-GATK变异检测工具包 | BWA-GATK Variant Calling Toolkit

功能: 完整的FASTQ到VCF流程，包括自动索引构建、GVCF模式变异检测和双重过滤
Features: Complete FASTQ-to-VCF pipeline with automatic index building, GVCF mode variant calling and dual filtering

作者 | Author: Claude
版本 | Version: v1.0
日期 | Date: 2025-10-06
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import BWAGATKPipeline
from .config import PipelineConfig

__all__ = ['BWAGATKPipeline', 'PipelineConfig']
