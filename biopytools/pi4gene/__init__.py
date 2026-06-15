"""
基因分组核苷酸多样性计算|Nucleotide Diversity Calculation per Gene Group
功能: 按分组提取序列、MAFFT比对、计算pi|
Features: Extract sequences by group, MAFFT alignment, calculate pi
作者|Author: Xiang LI
版本|Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import Pi4GeneAnalyzer
from .config import Pi4GeneConfig

__all__ = [
    'Pi4GeneAnalyzer',
    'Pi4GeneConfig'
]
