"""
基因组装配统计模块|Genome Assembly Statistics Module
"""

from .main import AssemblyStatsRunner
from .config import AssemblyStatsConfig
from .stats_analyzer import AssemblyStatsAnalyzer

__version__ = "1.0.0"

__all__ = [
    'AssemblyStatsRunner',
    'AssemblyStatsConfig',
    'AssemblyStatsAnalyzer'
]
