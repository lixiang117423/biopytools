"""
SignalP 6.0信号肽预测模块|SignalP 6.0 Signal Peptide Prediction Module
"""

from .main import SignalPPredictor, main
from .config import SignalPConfig

__version__ = "1.0.0"

__all__ = [
    'SignalPPredictor',
    'SignalPConfig',
    'main'
]
