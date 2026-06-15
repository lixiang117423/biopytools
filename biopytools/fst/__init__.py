"""
Fst计算工具包|Fst Calculation Toolkit
功能: 计算群体间遗传分化系数(Fst)，支持VCF输入和灵活的质量控制|
Features: Calculate genetic differentiation coefficient (Fst) between populations, supporting VCF input and flexible quality control
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-03-29

使用示例|Usage Examples:
    from biopytools.fst import FstAnalyzer, FstConfig

    # 创建分析器|Create analyzer
    analyzer = FstAnalyzer(
        vcf_file="variants.vcf",
        pop_file="population.txt",
        output_dir="./fst_output"
    )

    # 运行分析|Run analysis
    results = analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import FstAnalyzer
from .config import FstConfig
from .fst_calculator import FstCalculator
from .results_processor import FstResultsProcessor

__all__ = [
    'FstAnalyzer',
    'FstConfig',
    'FstCalculator',
    'FstResultsProcessor'
]
