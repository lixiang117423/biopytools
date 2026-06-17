"""
FASTA序列ID重命名模块|FASTA Sequence ID Renamer Module

用于简化NCBI下载的基因组序列ID，使其更适合后续分析使用
Simplify NCBI genome sequence IDs for better usability in downstream analysis
"""

from .main import main as fasta_id_renamer_main
from .config import FastaIDRenamerConfig
from .processor import FastaIDRenamer

__version__ = "1.0.0"

__all__ = [
    'fasta_id_renamer_main',
    'FastaIDRenamerConfig',
    'FastaIDRenamer'
]
