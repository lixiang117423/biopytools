"""
BUSCO质量评估分析工具包|BUSCO Quality Assessment Analysis Toolkit
功能|Features: 基因组/转录组/蛋白质序列的BUSCO质量评估完整流程，支持批处理和结果汇总|
Complete pipeline for BUSCO quality assessment of genome/transcriptome/protein sequences, supporting batch processing and result summarization
"""

__version__ = "1.0.0"

from .main import BUSCOAnalyzer
from .config import BUSCOConfig
from .analyzer import BUSCORunner
from .results import ResultsProcessor, SummaryGenerator
from .utils import BUSCOLogger, CommandRunner, FileManager

__all__ = [
    'BUSCOAnalyzer',
    'BUSCOConfig',
    'BUSCORunner',
    'ResultsProcessor',
    'SummaryGenerator',
    'BUSCOLogger',
    'CommandRunner',
    'FileManager'
]
