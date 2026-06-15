"""基于seqkit的FASTQ文件统计模块|FASTQ File Statistics Module Based on seqkit"""

from .main import FastqStatsRunner
from .config import FastqStatsConfig
from .calculator import FastqStatsCalculator

__version__ = "1.0.0"

__all__ = [
    'FastqStatsRunner',
    'FastqStatsConfig',
    'FastqStatsCalculator'
]
