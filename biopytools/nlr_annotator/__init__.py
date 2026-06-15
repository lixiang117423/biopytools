"""
NLR-Annotator模块|NLR-Annotator Module

运行NLR-Annotator从CDS序列中预测NLR基因
Run NLR-Annotator to predict NLR genes from CDS sequences
"""

from .config import NLRAnnotatorConfig
from .main import main

__version__ = "1.0.0"

__all__ = ['NLRAnnotatorConfig', 'main']
