"""
VCF SNP/INDEL过滤工具包 | VCF SNP/INDEL Filtering Toolkit
功能: VCF文件中SNP和INDEL的自动分离、过滤和统计 | 
Features: Automatic separation, filtering and statistics of SNPs and INDELs in VCF files
作者 | Author: Claude
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-10-11

使用示例 | Usage Examples:
    from biopytools.filter_snp_indel import VCFFilterAnalyzer, FilterConfig
    
    # 创建分析器 | Create analyzer
    analyzer = VCFFilterAnalyzer(
        vcf_file="variants.vcf",
        output_dir="filtered_output",
        threads=88
    )
    
    # 运行过滤 | Run filtering
    analyzer.run_filtering()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import VCFFilterAnalyzer
from .config import FilterConfig

__all__ = ['VCFFilterAnalyzer', 'FilterConfig']
