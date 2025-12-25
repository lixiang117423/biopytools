"""
GenomeThreader 基因预测工具包 | GenomeThreader Gene Prediction Toolkit
功能: 基于相似性比对的基因结构预测，支持cDNA/EST和蛋白质序列比对 | 
Features: Similarity-based gene structure prediction supporting cDNA/EST and protein sequence alignments
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-08-26

使用示例 | Usage Examples:
    from biopytools.genomethreader import GenomeThreaderAnalyzer, GTHConfig
    
    # 创建分析器 | Create analyzer
    analyzer = GenomeThreaderAnalyzer(
        genomic_file="genome.fa",
        cdna_file="cdna.fa",
        output_dir="gth_results"
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import GenomeThreaderAnalyzer
from .config import GTHConfig

__all__ = ['GenomeThreaderAnalyzer', 'GTHConfig']
