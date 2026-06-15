"""
重复序列屏蔽模块|Repeat Masking Module

该模块用于基因组重复序列的识别和屏蔽|This module is for genome repeat identification and masking

主要功能|Main features:
- RepeatModeler: 从头识别重复序列|De novo repeat identification
- RepeatMasker: 重复序列屏蔽|Repeat sequence masking
- 支持Dfam/Repbase数据库|Support for Dfam/Repbase databases
"""

from .main import main
from .config import RepeatMaskConfig
from .core import RepeatMaskPipeline, RepeatModelerRunner, RepeatMaskerRunner
from .utils import RepeatMaskLogger, CommandRunner, get_genome_stats

__version__ = "1.0.0"

__all__ = [
    'main',
    'RepeatMaskConfig',
    'RepeatMaskPipeline',
    'RepeatModelerRunner',
    'RepeatMaskerRunner',
    'RepeatMaskLogger',
    'CommandRunner',
    'get_genome_stats'
]
