"""
ANNOVAR基因注释工具包|ANNOVAR Gene Annotation Toolkit
功能: VCF变异注释的完整流程，支持GFF3转换和基因注释|
Features: Complete pipeline for VCF variant annotation, supporting GFF3 conversion and gene annotation
作者|Author: Xiang LI  
版本|Version: v10 - 模块化重构版|Modular refactored version
日期|Date: 2025-07-10

使用示例|Usage Examples:
    from biopytools.annovar import ANNOVARAnnotator, ANNOVARConfig

    # 创建注释器|Create annotator
    annotator = ANNOVARAnnotator(
        gff3_file="annotation.gff3",
        genome_file="genome.fa",
        vcf_file="variants.vcf",
        build_ver="OV",
        annovar_path="/path/to/annovar"
    )

    # 运行注释|Run annotation
    annotator.run_full_pipeline()

    # 处理注释结果（可选）| Process annotation results (optional)
    # annotator.process_annotation_results()
    #
    # 或者单独处理特定结果文件|Or process specific result files separately
    # annotator.process_exonic_results_only()
    # annotator.process_all_results_only()
"""

__version__ = "10.0.0"
__author__ = "Xiang LI"

from .main import ANNOVARAnnotator
from .config import ANNOVARConfig
from .results_processor import ANNOVARResultsProcessor, ExonicVariantProcessor, AllVariantProcessor

__all__ = ['ANNOVARAnnotator', 'ANNOVARConfig', 'ANNOVARResultsProcessor',
           'ExonicVariantProcessor', 'AllVariantProcessor']
