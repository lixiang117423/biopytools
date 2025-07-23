"""
VCF文件筛选工具包 | VCF File Filtering Toolkit
功能: VCF文件筛选的完整流程，支持多种过滤条件和格式转换 | 
Features: Complete pipeline for VCF file filtering, supporting various filtering conditions and format conversion
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-07-23

使用示例 | Usage Examples:
    from biopytools.vcf_filter import VCFFilterMain, VCFFilterConfig
    
    # 创建分析器 | Create analyzer
    analyzer = VCFFilterMain(
        vcf_file="variants.vcf",
        output_file="filtered.vcf",
        chr_name="1",
        start=1000,
        end=2000,
        convert_format=True,
        skip_validation=True  # 跳过验证以提高性能
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import VCFFilterMain, filter_vcf_file
from .config import VCFFilterConfig

__all__ = ['VCFFilterMain', 'VCFFilterConfig', 'filter_vcf_file']
