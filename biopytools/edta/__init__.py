"""
🌾 EDTA植物基因组TE注释工具包 | EDTA Plant Genome TE Annotation Toolkit
功能: 基于EDTA的植物基因组转座元件鉴定、分类和注释的完整流程
作者 | Author: Claude  
版本 | Version: v1.0
日期 | Date: 2025-08-26
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import EDTAAnalyzer
from .config import EDTAConfig

__all__ = ['EDTAAnalyzer', 'EDTAConfig']
