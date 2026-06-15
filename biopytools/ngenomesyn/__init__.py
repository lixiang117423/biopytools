"""
NGenomeSyn Module
NGenomeSyn基因组共线性分析模块
"""

from .config import NGenomeSynConfig
from .main import NGenomeSynAnalyzer, setup_logger
from .utils import NGenomeSynLogger

__all__ = [
    'NGenomeSynConfig',
    'NGenomeSynAnalyzer',
    'setup_logger',
    'NGenomeSynLogger'
]
