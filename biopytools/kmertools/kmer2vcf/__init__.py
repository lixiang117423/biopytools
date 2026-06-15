"""
Kmer丰度转VCF模块|Kmer Abundance to VCF Converter Module
"""

from .main import main as kmer2vcf_main
from .config import Kmer2VcfConfig
from .calculator import KmerToVcfConverter

__version__ = "1.0.0"

__all__ = [
    'kmer2vcf_main',
    'Kmer2VcfConfig',
    'KmerToVcfConverter'
]
