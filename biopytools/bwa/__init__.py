"""
BWA全基因组比对工具包|BWA Whole Genome Alignment Toolkit
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import BWAAligner
from .config import AlignConfig

__all__ = ['BWAAligner', 'AlignConfig']
