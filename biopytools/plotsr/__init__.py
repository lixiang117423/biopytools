"""
多基因组共线性可视化工具|Multi-Genome Synteny Visualization Tool

此模块封装了minimap2 + SyRI + PlotSR完整流程，实现一键式多基因组共线性分析
This module wraps minimap2 + SyRI + PlotSR pipeline for one-click multi-genome synteny analysis
"""

from .main import main as plotsr_main
from .config import PlotSRConfig
from .pipeline import PlotSRPipeline

__version__ = "1.0.0"

__all__ = [
    'plotsr_main',
    'PlotSRConfig',
    'PlotSRPipeline'
]
