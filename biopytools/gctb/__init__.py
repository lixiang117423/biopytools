"""
GCTB全基因组复杂性状贝叶斯分析模块|GCTB Genome-wide Complex Trait Bayesian Analysis Module
"""

from .main import GCTBRunner
from .config import GCTBConfig

__version__ = "1.0.0"

__all__ = [
    'GCTBRunner',
    'GCTBConfig'
]
