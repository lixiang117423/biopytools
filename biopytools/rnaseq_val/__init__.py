"""
转录组验证注释|Transcriptome Validation for Genome Annotation

用二代/三代转录组数据验证基因组注释转录本结构，自动分级校正
"""

from .main import RnaseqValPipeline
from .config import RnaseqValConfig

__version__ = "1.0.0"

__all__ = [
    'RnaseqValPipeline',
    'RnaseqValConfig',
]
