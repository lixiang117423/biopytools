"""
PopLDdecay封装模块|PopLDdecay Wrapper Module

用于计算和可视化连锁不平衡衰减|Calculate and visualize linkage disequilibrium decay
"""

from .main import PopLDdecayRunner
from .ld_threshold import LDThresholdRecommender, recommend_ld_threshold

__all__ = [
    'PopLDdecayRunner',
    'LDThresholdRecommender',
    'recommend_ld_threshold'
]
__version__ = '1.1.0'
