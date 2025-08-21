"""
🧬 BLAST比对分析工具包 | BLAST Alignment Analysis Toolkit
功能: 多序列文件与目标序列的BLAST比对分析完整流程
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本
日期 | Date: 2025-08-21
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import BLASTAnalyzer
from .config import BLASTConfig

__all__ = ['BLASTAnalyzer', 'BLASTConfig']
