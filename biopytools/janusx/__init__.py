"""
JanusX GWAS和基因组选择分析工具包|JanusX GWAS and Genomic Selection Analysis Toolkit

JanusX是一个高性能的全基因组关联分析(GWAS)和基因组选择(GS)工具包
JanusX is a high-performance Genome-Wide Association Study (GWAS) and Genomic Selection (GS) toolkit
"""

from .main import main
from .config import JanusXGWASConfig, JanusXGSConfig, JanusXPCAConfig, JanusXPostGWASConfig

__version__ = "1.0.0"

__all__ = [
    'main',
    'JanusXGWASConfig',
    'JanusXGSConfig',
    'JanusXPCAConfig',
    'JanusXPostGWASConfig'
]
