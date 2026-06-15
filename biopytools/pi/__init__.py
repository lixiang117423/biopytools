"""
核苷酸多样性计算工具包|Nucleotide Diversity Calculation Toolkit
功能: 使用vcftools计算群体内核苷酸多样性(pi)，输出汇总表|
Features: Calculate within-population nucleotide diversity (pi) using vcftools, output summary table
作者|Author: Xiang LI
版本|Version: 2.0.0
日期|Date: 2026-04-13

使用示例|Usage Examples:
    from biopytools.pi import PiAnalyzer, PiConfig

    # 创建分析器|Create analyzer
    analyzer = PiAnalyzer(
        vcf_file="variants.vcf.gz",
        pop_file="population.txt",
        genome="reference.fasta",
        output_dir="./pi_output"
    )

    # 运行分析|Run analysis
    analyzer.run()
"""

__version__ = "2.0.0"
__author__ = "Xiang LI"

from .main import PiAnalyzer
from .config import PiConfig

__all__ = [
    'PiAnalyzer',
    'PiConfig'
]
