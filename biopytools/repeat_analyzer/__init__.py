"""
重复序列分析工具包 | Repeat Sequence Analysis Toolkit
功能: 植物基因组重复序列和转座元件的完整分析流程 | 
Features: Complete pipeline for plant genome repeat sequences and transposable elements analysis
作者 | Author: Claude
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-08-26

🚀 使用示例 | Usage Examples:
    from biopytools.repeat_analyzer import RepeatAnalyzer
    
    # 创建分析器 | Create analyzer
    analyzer = RepeatAnalyzer(
        genome_file="genome.fasta",
        output_dir="repeat_results",
        threads=88
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import RepeatAnalyzer
from .config import RepeatConfig

__all__ = ['RepeatAnalyzer', 'RepeatConfig']
