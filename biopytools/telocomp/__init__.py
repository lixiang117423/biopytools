"""
TeloComp端粒鉴定和可视化工具包|TeloComp Telomere Identification and Visualization Toolkit
功能: 端粒序列鉴定、提取、补全和可视化分析|
Features: Telomere sequence identification, extraction, complementation and visualization
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-01-15

使用示例|Usage Examples:
    from biopytools.telocomp import TeloCompAnalyzer, TeloCompConfig

    # 创建分析器|Create analyzer
    analyzer = TeloCompAnalyzer(
        genome="genome.fa",
        ont="ont.fastq.gz",
        hifi="hifi.fastq.gz",
        output_dir="telomere_results"
    )

    # 运行端粒鉴定流程|Run telomere identification pipeline
    analyzer.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import TeloCompAnalyzer
from .config import TeloCompConfig

__all__ = ['TeloCompAnalyzer', 'TeloCompConfig']
