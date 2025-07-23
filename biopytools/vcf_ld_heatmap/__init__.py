"""
VCF LD热图生成工具包 | VCF LD Heatmap Generation Toolkit
功能: VCF文件连锁不平衡热图生成的完整流程，支持多种过滤和可视化选项 | 
Features: Complete pipeline for VCF linkage disequilibrium heatmap generation, supporting various filtering and visualization options
作者 | Author: Xiang LI  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-07-23

使用示例 | Usage Examples:
    from biopytools.vcf_ld_heatmap import LDHeatmapAnalyzer, LDHeatmapConfig
    
    # 创建分析器 | Create analyzer
    analyzer = LDHeatmapAnalyzer(
        vcf_file="variants.vcf",
        output_file="ld_heatmap.png",
        maf=0.05,
        max_snps=500
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import LDHeatmapAnalyzer
from .config import LDHeatmapConfig

__all__ = ['LDHeatmapAnalyzer', 'LDHeatmapConfig']
