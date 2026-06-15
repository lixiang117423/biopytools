"""
PanDepth覆盖度计算工具包|PanDepth Coverage Calculation Toolkit
功能: 基于PanDepth的超快速基因组覆盖度计算|Features: Ultra-fast genome coverage calculation based on PanDepth
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-06

使用示例|Usage Examples:
    from biopytools.pandepth import PanDepthCalculator, PanDepthConfig

    # 创建计算器|Create calculator
    calculator = PanDepthCalculator(
        input="sample.bam",
        output_dir="coverage_results",
        threads=12
    )

    # 运行覆盖度计算|Run coverage calculation
    calculator.run()

    # 基因覆盖度计算|Gene coverage calculation
    calculator = PanDepthCalculator(
        input="sample.bam",
        gff_file="genes.gff",
        output_dir="gene_coverage",
        threads=12
    )
    calculator.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import PanDepthCalculator
from .config import PanDepthConfig
from .results import PanDepthResultsMerger

__all__ = ['PanDepthCalculator', 'PanDepthConfig', 'PanDepthResultsMerger']
