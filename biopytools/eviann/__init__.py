"""
EviAnn基因组注释模块|EviAnn Genome Annotation Module
"""

from .main import EviAnnotator
from .config import EviAnnConfig

__version__ = "1.0.0"

__all__ = [
    'EviAnnotator',
    'EviAnnConfig'
]
