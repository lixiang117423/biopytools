"""
🧬 K-mer丰度分析包初始化 | K-mer Abundance Analysis Package Initialization
"""

__version__ = "1.0.0"
__author__ = "biopytools team"

from .main import KmerCountAnalyzer, main
from .config import KmerCountConfig

__all__ = ['KmerCountAnalyzer', 'KmerCountConfig', 'main']
