"""
VCF系统发育分析工具包 | VCF Phylogenetic Analysis Toolkit
功能: VCF文件系统发育分析的完整流程，支持距离矩阵计算和NJ树构建 | 
Features: Complete pipeline for VCF phylogenetic analysis, supporting distance matrix calculation and NJ tree construction
作者 | Author: Claude  
版本 | Version: v2.0 - 简化scikit-bio版本 | Simplified scikit-bio version
日期 | Date: 2025-07-25

使用示例 | Usage Examples:
    from biopytools.vcf_phylo import VCFPhyloAnalyzer, PhyloConfig
    
    # 创建分析器 | Create analyzer
    analyzer = VCFPhyloAnalyzer(
        vcf_file="variants.vcf",
        output_prefix="phylo_results"
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "2.0.0"
__author__ = "Claude"

from .main import VCFPhyloAnalyzer
from .config import PhyloConfig

__all__ = ['VCFPhyloAnalyzer', 'PhyloConfig']
