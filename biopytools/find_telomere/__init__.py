"""
端粒识别分析模块|Telomere Identification Module

该模块提供了端粒重复序列识别和分析的功能，基于 tidk (Telomere Identification toolKit) 工具。
支持探索、查找、搜索和可视化端粒重复序列。

This module provides telomeric repeat identification and analysis functionality based on the
tidk (Telomere Identification toolKit) tool. Supports exploration, finding, searching,
and visualization of telomeric repeats.
"""

from .main import TelomereFinder

__version__ = '1.0.0'
__all__ = ['TelomereFinder']
