"""基因密度计算|Gene Density Calculation

按固定大小窗口统计每条染色体各区间的基因数量与基因密度(基因/Mb)
|Count genes and gene density (genes/Mb) per fixed-size window along each chromosome
"""

from .config import GeneDensityConfig
from .calculator import GeneDensityCalculator

__version__ = "1.0.0"

__all__ = [
    'GeneDensityConfig',
    'GeneDensityCalculator',
]
