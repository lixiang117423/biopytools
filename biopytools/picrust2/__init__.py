"""
PICRUSt2微生物群落功能丰度预测工具|PICRUSt2 Microbial Community Functional Abundance Prediction
功能: 基于16S rRNA标记基因序列预测微生物群落功能丰度和代谢通路|
Features: Predict functional abundances and metabolic pathways from marker gene sequences
"""

__version__ = "1.0.0"

from .pipeline import Picrust2Pipeline
from .config import Picrust2Config

__all__ = ['Picrust2Pipeline', 'Picrust2Config']
