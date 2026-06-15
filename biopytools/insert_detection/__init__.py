"""插入检测模块|Insert detection module"""

from .main import main
from .config import InsertDetectionConfig
from .detector import InsertDetector
from .output import write_results, write_summary

__version__ = "1.0.0"

__all__ = [
    'main',
    'InsertDetectionConfig',
    'InsertDetector',
    'write_results',
    'write_summary'
]
