"""
🧬 基因组共线性可视化工具包 | Genome Synteny Visualization Toolkit
"""

__version__ = "1.0.0"
__author__ = "🧑‍💻 Xiang LI"

from .main import GenomeSynAnalyzer, main
from .config import GenomeSynConfig

__all__ = ['GenomeSynAnalyzer', 'GenomeSynConfig', 'main']
