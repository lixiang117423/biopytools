"""
Minibwa短读长比对工具包|Minibwa Short-read Alignment Toolkit

基于lh3的minibwa (https://github.com/lh3/minibwa) 封装，提供短读长/Hi-C/BS-seq模式比对
|Wrapper for lh3's minibwa, supporting short-read/Hi-C/BS-seq alignment modes
"""

from .main import MinibwaRunner
from .config import MinibwaConfig

__version__ = "1.0.0"

__all__ = [
    'MinibwaRunner',
    'MinibwaConfig',
]
