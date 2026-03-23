"""
VCF PCA分析工具包|VCF PCA Analysis Toolkit
功能: VCF文件主成分分析的完整流程，支持质控、可视化和样本信息整合|
Features: Complete pipeline for VCF principal component analysis, supporting QC, visualization and sample info integration
作者|Author: Xiang LI  
版本|Version: v1.0 - 模块化版本|Modular version
日期|Date: 2025-07-22

使用示例|Usage Examples:
    from biopytools.vcf_pca import VCFPCAAnalyzer, PCAConfig
    
    # 创建分析器|Create analyzer
    analyzer = VCFPCAAnalyzer(
        vcf_file="variants.vcf",
        output_dir="pca_results",
        components=10,
        plot=True
    )
    
    # 运行分析|Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import VCFPCAAnalyzer
from .config import PCAConfig

__all__ = ['VCFPCAAnalyzer', 'PCAConfig']
