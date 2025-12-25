"""
Filter ANNOVAR 工具包 | Filter ANNOVAR Toolkit
基因区域变异提取的完整流程 | Complete pipeline for gene region variant extraction
"""

__version__ = "1.0.0"
__author__ = "Bioinformatics Team"

from .main import FilterAnnovarAnalyzer
from .config import FilterConfig

__all__ = ['FilterAnnovarAnalyzer', 'FilterConfig']
