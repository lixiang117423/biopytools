"""
VCF2Gene变异注释工具包|VCF2Gene Variant Annotation Toolkit
功能: 基于GFF注释文件对VCF变异进行基因区域注释|
Features: Annotate VCF variants with gene regions based on GFF annotation
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-05

使用示例|Usage Examples:
    from biopytools.vcf2gene import VCF2GeneRunner, VCF2GeneConfig

    # 创建注释器|Create annotator
    runner = VCF2GeneRunner(
        vcf_file="variants.vcf",
        gff_file="annotation.gff",
        output_file="annotated_variants.txt"
    )

    # 运行注释|Run annotation
    runner.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import VCF2GeneRunner
from .config import VCF2GeneConfig
from .calculator import VariantAnnotator, GFFParser

__all__ = [
    'VCF2GeneRunner',
    'VCF2GeneConfig',
    'VariantAnnotator',
    'GFFParser'
]
