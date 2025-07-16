"""
VCF基因型提取工具包 | VCF Genotype Extraction Toolkit
功能: 从VCF文件中提取基因型信息，支持多种输出格式和过滤选项 | 
Features: Extract genotype information from VCF files, supporting multiple output formats and filtering options
作者 | Author: Claude  
版本 | Version: v1.0.0
日期 | Date: 2025-07-16

使用示例 | Usage Examples:
    from biopytools.vcf_genotype_extractor import VCFGenotypeExtractor, VCFConfig
    
    # 创建提取器 | Create extractor
    extractor = VCFGenotypeExtractor(
        vcf_file="variants.vcf.gz",
        output_prefix="genotype_results",
        samples="all",
        output_type="txt"
    )
    
    # 运行提取 | Run extraction
    extractor.run_extraction()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import VCFGenotypeExtractor
from .config import VCFConfig

__all__ = ['VCFGenotypeExtractor', 'VCFConfig']
