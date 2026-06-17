"""
HiC-Pro质量控制评估模块|HiC-Pro Quality Control Assessment Module

该模块用于解析HiC-Pro输出并评估Hi-C数据质量
This module parses HiC-Pro output and assesses Hi-C data quality
"""

from .config import HiCProQCConfig
from .calculator import HiCProQCCalculator
from .main import HiCProQCRunner, main as hicpro_qc_main

__all__ = [
    'HiCProQCConfig',
    'HiCProQCCalculator',
    'HiCProQCRunner',
    'hicpro_qc_main'
]
