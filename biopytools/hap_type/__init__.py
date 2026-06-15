"""
单倍型分析模块|Haplotype Analysis Module

从VCF文件提取基因单倍型，支持geneHapR兼容输出|Extract gene haplotypes from VCF, geneHapR-compatible output
"""

from .config import HapTypeConfig
from .utils import HapTypeLogger
from .builder import HaplotypeBuilder, AlleleEncoder
from .main import main

__all__ = [
    'HapTypeConfig',
    'HapTypeLogger',
    'HaplotypeBuilder',
    'AlleleEncoder',
    'main'
]
