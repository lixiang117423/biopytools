"""
基因组挂载率统计工具包|Genome Mount Rate Statistics Toolkit
功能: 计算FASTA文件中前N条或最长N条序列占总基因组长度的百分比|
Features: Calculate the percentage of top N or longest N sequences in total genome length
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-10

使用示例|Usage Examples:
    from biopytools.genome_mount_rate import GenomeMountRateCalculator, GenomeMountRateConfig

    # 创建计算器|Create calculator
    calculator = GenomeMountRateCalculator(
        fasta_file="genome.fa",
        number=10,
        sort_by_length=True
    )

    # 运行计算|Run calculation
    results = calculator.calculate()
    print(f"挂载率|Mount rate: {results['percentage']:.2f}%")
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import GenomeMountRateRunner
from .config import GenomeMountRateConfig
from .calculator import GenomeMountRateCalculator

__all__ = ['GenomeMountRateRunner', 'GenomeMountRateConfig', 'GenomeMountRateCalculator']
