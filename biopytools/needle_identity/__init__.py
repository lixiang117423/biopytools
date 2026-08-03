"""序列两两identity计算模块(EMBOSS needle)|Pairwise Sequence Identity Module (EMBOSS needle)"""

from .calculator import NeedleIdentityCalculator
from .config import NeedleIdentityConfig
from .utils import NeedleIdentityLogger

__version__ = "1.0.0"

__all__ = [
    'NeedleIdentityConfig',
    'NeedleIdentityCalculator',
    'NeedleIdentityLogger',
]
