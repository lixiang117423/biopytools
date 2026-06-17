"""
TGS-GapCloser Gap填充工具包|TGS-GapCloser Gap Filling Toolkit
功能: 使用三代测序数据填充基因组组装中的Gap|Features: Fill gaps in genome assembly using TGS long reads
版本|Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import TGSGapCloser
from .config import TGSGapCloserConfig

__all__ = ['TGSGapCloser', 'TGSGapCloserConfig']
