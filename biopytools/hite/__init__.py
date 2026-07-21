"""
HiTE转座子检测与注释模块|HiTE Transposon Detection and Annotation Module

包含 HiTE 单基因组分析(singularity 直接挂载)和 panHiTE 群体基因组分析(保留)
Includes HiTE single-genome analysis (singularity direct-mount) and panHiTE (preserved)
"""

__version__ = "1.1.0"

from .config import HiteConfig, PanHiteConfig
from .hite_runner import HiteRunner
from .results import HiteResultsProcessor, PanHiteResultsProcessor
from .panhite_analyzer import PanHiteAnalyzer

__all__ = [
    'HiteConfig',
    'HiteRunner',
    'HiteResultsProcessor',
    'PanHiteConfig',
    'PanHiteResultsProcessor',
    'PanHiteAnalyzer',
]
