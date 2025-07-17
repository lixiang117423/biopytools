"""
VCF基因型统计工具包 | VCF Genotype Statistics Toolkit
功能: VCF文件中每个样品的基因型统计分析，包括杂合率和纯合率计算 | 
Features: Genotype statistics analysis for each sample in VCF files, including heterozygosity and homozygosity rates

支持的基因型格式 | Supported genotype formats:
- 未定相基因型: 0/0, 0/1, 1/1, ./. | Unphased genotypes: 0/0, 0/1, 1/1, ./.
- 已定相基因型: 0|0, 0|1, 1|1, .|. | Phased genotypes: 0|0, 0|1, 1|1, .|.
- 多等位基因: 0/2, 1/2, 2/2等 | Multi-allelic: 0/2, 1/2, 2/2, etc.

作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-07-17

使用示例 | Usage Examples:
    from biopytools.vcf_stats import VCFStatsAnalyzer, VCFStatsConfig
    
    # 创建分析器 | Create analyzer
    analyzer = VCFStatsAnalyzer(
        vcf_file="variants.vcf",
        output_dir="vcf_stats_output"
    )
    
    # 运行统计分析 | Run statistics analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import VCFStatsAnalyzer
from .config import VCFStatsConfig

__all__ = ['VCFStatsAnalyzer', 'VCFStatsConfig']
