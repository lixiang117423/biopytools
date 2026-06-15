"""
Protein Stats蛋白质理化性质分析工具包|Protein Stats Protein Properties Analysis Toolkit
功能: 计算蛋白质序列的理化性质（长度、分子量、等电点等）|
Features: Calculate protein sequence properties (length, MW, pI, etc.)
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-28

使用示例|Usage Examples:
    from biopytools.protein_stats import ProteinStatsPipeline

    # 创建分析器|Create pipeline
    pipeline = ProteinStatsPipeline(
        protein_fasta="proteins.fa",
        output_file="stats.tsv"
    )

    # 运行分析|Run analysis
    pipeline.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import ProteinStatsPipeline
from .config import ProteinStatsConfig
from .analyzer import ProteinStatsAnalyzer

__all__ = [
    'ProteinStatsPipeline',
    'ProteinStatsConfig',
    'ProteinStatsAnalyzer'
]
