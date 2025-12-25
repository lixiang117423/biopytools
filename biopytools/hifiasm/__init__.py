"""
HiFiasm基因组组装分析模块 | HiFiasm Genome Assembly Analysis Module
"""

__version__ = "1.0.0"
__author__ = "Bioinformatics Team"
__description__ = "HiFiasm基因组组装完整流水线 | HiFiasm Genome Assembly Complete Pipeline"

from .config import HifiasmConfig
from .main import HifiasmAnalyzer
from .utils import HifiasmLogger, CommandRunner, check_dependencies
from .assembly import HifiasmAssembler, GFAConverter
from .quality_assessment import BUSCOAssessor, QUASTAssessor, HaplotypeAnalyzer
from .data_processing import StatisticsCalculator, FormatConverter

__all__ = [
    'HifiasmConfig',
    'HifiasmAnalyzer', 
    'HifiasmLogger',
    'CommandRunner',
    'check_dependencies',
    'HifiasmAssembler',
    'GFAConverter',
    'BUSCOAssessor',
    'QUASTAssessor',
    'HaplotypeAnalyzer',
    'StatisticsCalculator',
    'FormatConverter'
]