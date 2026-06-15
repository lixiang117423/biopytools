"""Hi-C热图分析模块|Hi-C heatmap analysis module (HiCPro + PlotHiC)"""

from .main import main
from .config import HiCProConfig
from .hicpro_pipeline import HiCProPipeline
from .utils import HiCLogger

__version__ = "3.0.0"

__all__ = [
    'main',
    'HiCProConfig',
    'HiCProPipeline',
    'HiCLogger'
]
