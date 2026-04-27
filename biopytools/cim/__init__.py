"""
CIM分析模块|R/qtl Composite Interval Mapping (CIM) Analysis Module
"""

from .main import main as cim_main
from .config import CIMConfig

__version__ = "1.0.0"

__all__ = [
    'cim_main',
    'CIMConfig'
]
