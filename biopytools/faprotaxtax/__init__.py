"""
FAPROTAX微生物群落功能注释工具|FAPROTAX Functional Annotation of Prokaryotic Taxa
功能: 基于分类学注释将微生物群落OTU/ASV表转换为功能丰度表|
Features: Convert taxonomic microbial community profiles to functional abundance profiles
"""

from .pipeline import FaprotaxtaxPipeline
from .config import FaprotaxtaxConfig

__version__ = "1.0.0"

__all__ = [
    'FaprotaxtaxPipeline',
    'FaprotaxtaxConfig'
]
