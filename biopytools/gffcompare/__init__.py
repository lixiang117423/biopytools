"""
GFF/GTF文件两两比较分析模块|GFF/GTF Pairwise Comparison Analysis Module
功能|Features: 使用gffcompare对多个GFF/GTF文件进行两两双向比较|
Pairwise bidirectional comparison of multiple GFF/GTF files using gffcompare
"""

__version__ = "1.0.0"

from .main import GffComparePairwise
from .config import GffCompareConfig
from .utils import GffCompareLogger, CommandRunner

__all__ = [
    'GffComparePairwise',
    'GffCompareConfig',
    'GffCompareLogger',
    'CommandRunner'
]
