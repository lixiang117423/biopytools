"""INDEL分子标记开发模块|INDEL Marker Module"""

from .main import IndelMarkerRunner
from .config import IndelMarkerConfig
from .samplesheet import SampleInfo, SamplesheetParser

__version__ = "1.0.0"

__all__ = [
    'IndelMarkerRunner',
    'IndelMarkerConfig',
    'SampleInfo',
    'SamplesheetParser',
]
