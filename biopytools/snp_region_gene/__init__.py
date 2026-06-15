"""SNP区域基因提取模块|SNP Region Gene Extractor Module"""

from .main import main as snp_region_gene_main
from .config import SnpRegionConfig
from .processor import SnpRegionProcessor

__version__ = "1.0.0"

__all__ = [
    'snp_region_gene_main',
    'SnpRegionConfig',
    'SnpRegionProcessor'
]
