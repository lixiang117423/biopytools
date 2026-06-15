"""
HiTE转座子检测与注释模块|HiTE Transposon Detection and Annotation Module

包含 HiTE 单基因组分析和 panHiTE 群体基因组分析功能
Includes HiTE single-genome analysis and panHiTE pan-genome analysis

HiTE: a fast and accurate dynamic boundary adjustment approach for
full-length Transposable Elements detection and annotation

Author: BioPyTools
Version: 1.0.0
"""

__version__ = "1.0.0"

from .hite_analyzer import HiteAnalyzer
from .panhite_analyzer import PanHiteAnalyzer
from .config import HiteConfig, PanHiteConfig
from .results import HiteResultsProcessor, PanHiteResultsProcessor

__all__ = [
    'HiteAnalyzer',
    'PanHiteAnalyzer',
    'HiteConfig',
    'PanHiteConfig',
    'HiteResultsProcessor',
    'PanHiteResultsProcessor'
]
