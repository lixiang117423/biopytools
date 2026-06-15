"""
OcBSA - BSA分析工具套件|OcBSA - BSA Analysis Tool Suite

F1群体DHHP分析、F2/RILs群体SNP-index/ED分析、BSA结果可视化、候选区域引物设计
"""

__version__ = "1.0.0"
__all__ = ['OcbsaCalculator', 'OcbsaConfig']


def __getattr__(name):
    """延迟导入重度依赖模块|Lazy import modules with heavy dependencies"""
    if name == 'OcbsaCalculator':
        from .core import OcbsaCalculator
        return OcbsaCalculator
    if name == 'OcbsaConfig':
        from .config import OcbsaConfig
        return OcbsaConfig
    raise AttributeError(f"module 'biopytools.ocbsa' has no attribute {name!r}")
