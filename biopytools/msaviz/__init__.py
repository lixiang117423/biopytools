"""
MSA可视化模块|MSA Visualization Module

多序列比对(MSA)可视化工具，基于matplotlib实现|Multiple Sequence Alignment (MSA) visualization tool implemented based on matplotlib
"""

from .msaviz import MsaViz
from .logger import MsaVizLogger

__version__ = "0.5.0"

__all__ = [
    'MsaViz',
    'MsaVizLogger',
]
