"""
EviAnn基因组注释模块|EviAnn Genome Annotation Module
"""

from .calculator import EviAnnRunner
from .config import EviAnnConfig

__version__ = "2.0.0"

__all__ = [
    'EviAnnRunner',
    'EviAnnConfig'
]
