"""
基因组重复序列分析工具包 | Genome Repeat Sequence Analysis Toolkit
功能: RepeatMasker、RepeatModeler、TRF的完整流程，支持重复序列鉴定和屏蔽 | 
Features: Complete pipeline for RepeatMasker, RepeatModeler, TRF, supporting repeat sequence identification and masking
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化重构版 | Modular refactored version
日期 | Date: 2025-07-20

使用示例 | Usage Examples:
    from biopytools.repeat_masker import RepeatAnalyzer, RepeatConfig
    
    # 创建分析器 | Create analyzer
    analyzer = RepeatAnalyzer(
        genome_file="genome.fasta",
        species="human",
        output_dir="repeat_analysis"
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import RepeatAnalyzer
from .config import RepeatConfig

__all__ = ['RepeatAnalyzer', 'RepeatConfig']
