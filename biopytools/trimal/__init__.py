"""trimal 多序列比对修剪工具包|trimal MSA Trimming Toolkit
功能: 封装 trimAl,对多序列比对做自动/手动修剪
|Features: Wrap trimAl for automated/manual trimming of multiple sequence alignments
"""

__version__ = "1.0.0"

from .main import TrimalRunner
from .config import TrimalConfig

__all__ = [
    'TrimalRunner',
    'TrimalConfig',
]
