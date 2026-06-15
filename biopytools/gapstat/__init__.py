"""
基因组Gap统计工具包|Genome Gap Statistics Toolkit
功能: 统计基因组FASTA文件中Gap的位置和长度|Features: Statistics gap positions and lengths in genome FASTA files
版本|Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import GapStat
from .config import GapStatConfig

__all__ = ['GapStat', 'GapStatConfig']
