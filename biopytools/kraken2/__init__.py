"""
Kraken2宏基因组分类工具包|Kraken2 Metagenomic Classification Toolkit
功能: 基于Kraken2的序列分类与物种丰度分析|
Features: Sequence classification and species abundance analysis based on Kraken2
作者|Author: Xiang LI
版本|Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import Kraken2Pipeline
from .config import Kraken2Config
from .analyzer import Kraken2Analyzer

__all__ = ['Kraken2Pipeline', 'Kraken2Config', 'Kraken2Analyzer']
