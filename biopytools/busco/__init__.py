"""
BUSCO质量评估分析工具包|BUSCO Quality Assessment Analysis Toolkit
功能|Features: 基因组/转录组/蛋白质序列的BUSCO质量评估完整流程，支持批处理和结果汇总|
Complete pipeline for BUSCO quality assessment of genome/transcriptome/protein sequences, supporting batch processing and result summarization
作者|Author: Claude
版本|Version: v1.0 - 模块化版本|Modular version
日期|Date: 2025-09-20

使用示例|Usage Examples:
    from biopytools.busco import BUSCOAnalyzer, BUSCOConfig

    # 创建分析器|Create analyzer
    analyzer = BUSCOAnalyzer(
        input_path="sequences.fa",
        lineage="brassicales_odb12",
        output_dir="busco_results",
        mode="genome",
        cpu=88
    )

    # 运行分析|Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import BUSCOAnalyzer
from .config import BUSCOConfig

__all__ = ['BUSCOAnalyzer', 'BUSCOConfig']
