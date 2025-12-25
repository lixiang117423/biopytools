"""
INDEL PAV分析工具包 | INDEL PAV Analysis Toolkit
功能: VCF文件INDEL存在缺失变异(PAV)分析的完整流程，支持多线程处理和结果可视化 | 
Features: Complete pipeline for INDEL Presence-Absence Variation (PAV) analysis from VCF files, supporting multi-threading and result visualization
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-01-15

使用示例 | Usage Examples:
    from biopytools.indel_pav import IndelPAVAnalyzer, PAVConfig
    
    # 创建分析器 | Create analyzer
    analyzer = IndelPAVAnalyzer(
        vcf_file="variants.vcf",
        output_file="indel_pav.txt",
        threads=88,
        min_length=5
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import IndelPAVAnalyzer
from .config import PAVConfig

__all__ = ['IndelPAVAnalyzer', 'PAVConfig']
