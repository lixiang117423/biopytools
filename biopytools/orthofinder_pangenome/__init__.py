"""
OrthoFinder泛基因组分析工具包 | OrthoFinder Pangenome Analysis Toolkit
功能: 基于OrthoFinder的泛基因组分析完整流程，支持核心基因组、附属基因组和特异基因组分类 | 
Features: Complete pangenome analysis pipeline based on OrthoFinder, supporting core, accessory and specific genome classification
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-08-28

使用示例 | Usage Examples:
    from biopytools.orthofinder_pangenome import OrthoFinderPangenomeAnalyzer, PangenomeConfig
    
    # 创建分析器 | Create analyzer
    analyzer = OrthoFinderPangenomeAnalyzer(
        input_dir="protein_sequences/",
        output_dir="pangenome_results",
        soft_threshold="all-1",
        search_program="blastp"
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import OrthoFinderPangenomeAnalyzer
from .config import PangenomeConfig

__all__ = ['OrthoFinderPangenomeAnalyzer', 'PangenomeConfig']
