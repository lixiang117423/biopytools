"""
Bismark甲基化分析流程工具包 | Bismark Methylation Pipeline Toolkit
"""

__version__ = "2.0.0"
__author__ = "Claude"

from .main import BismarkAnalyzer
from .config import BismarkConfig

__all__ = ['BismarkAnalyzer', 'BismarkConfig']
