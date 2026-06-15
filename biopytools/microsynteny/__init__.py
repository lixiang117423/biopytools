"""
微观共线性分析模块|Microsynteny Analysis Module

基于JCVI的自动化微观共线性分析和可视化工具
Automated microsynteny analysis and visualization tool based on JCVI
"""

from .main import main as microsynteny_main
from .config import MicrosyntenyConfig
from .analyzer import MicrosyntenyAnalyzer

__version__ = "1.0.0"

__all__ = [
    'microsynteny_main',
    'MicrosyntenyConfig',
    'MicrosyntenyAnalyzer'
]
