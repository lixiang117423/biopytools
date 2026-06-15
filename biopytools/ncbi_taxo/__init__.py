"""NCBI分类学注释模块|NCBI Taxonomy Annotation Module"""

from .main import NCBITaxoAnnotator
from .config import NCBITaxoConfig
from .processor import BlastTaxonomyProcessor
from .stats import TaxonomyStatsCalculator

__version__ = "1.0.0"

__all__ = [
    'NCBITaxoAnnotator',
    'NCBITaxoConfig',
    'BlastTaxonomyProcessor',
    'TaxonomyStatsCalculator'
]
