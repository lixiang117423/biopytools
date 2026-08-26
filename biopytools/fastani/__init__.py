"""fastANI 全基因组ANI计算工具包|fastANI whole-genome ANI toolkit
功能: 封装 FastANI,计算基因组间平均核苷酸一致性(ANI),输出矩阵与最近邻表
|Features: Wrap FastANI for whole-genome ANI, producing matrix and nearest-neighbor table
"""

__version__ = "1.0.0"

from .config import FastaniConfig
from .calculator import FastaniCalculator

__all__ = ['FastaniConfig', 'FastaniCalculator']
