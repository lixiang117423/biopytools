"""
甲基化分析工具包 | Methylation Analysis Pipeline Toolkit
"""

__version__ = "2.2.0"
__author__ = "Claude"

from .config import MethylationConfig
from .main import MethylationAnalyzer

__all__ = ["MethylationAnalyzer", "MethylationConfig"]
