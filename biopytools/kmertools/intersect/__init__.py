"""
Kmer交集提取模块|Kmer Intersection Module

从kmer矩阵中提取目标kmer的丰度信息|Extract kmer abundance from matrix based on target kmer list
支持正向和反向互补查询|Support forward and reverse complement queries
"""

from .main import main as intersect_main

__version__ = "1.0.0"

__all__ = ['intersect_main']
