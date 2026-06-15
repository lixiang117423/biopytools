"""
VCF SNP/INDEL过滤工具包|VCF SNP/INDEL Filtering Toolkit

功能: VCF文件中SNP和INDEL的自动分离、过滤和统计
Features: Automatic separation, filtering and statistics of SNPs and INDELs in VCF files

作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2024-12-30

使用示例|Usage Examples:
    from biopytools.filter_snp_indel import VCFFilterAnalyzer, FilterConfig

    # 创建分析器|Create analyzer
    analyzer = VCFFilterAnalyzer(
        vcf_file="variants.vcf",
        output_dir="filtered_output",
        threads=12
    )

    # 运行过滤|Run filtering
    analyzer.run_filtering()

命令行使用|Command Line Usage:
    python -m filter_snp_indel.main -i variants.vcf -o filtered_output/ -t 12
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import VCFFilterAnalyzer
from .config import FilterConfig

__all__ = ['VCFFilterAnalyzer', 'FilterConfig']
