"""区域单倍型分析(EasyHap 封装)|Regional haplotype analysis (EasyHap wrapper)"""

from .calculator import EasyHapCalculator
from .config import EasyHapConfig

__version__ = "1.0.0"

# 勿在此导出 main 函数(与子模块名冲突)|Do NOT export main here
__all__ = ["EasyHapConfig", "EasyHapCalculator", "__version__"]
