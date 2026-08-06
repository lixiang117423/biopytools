"""
FOF文件生成模块|FOF File Generation Module

生成样品名到文件路径的映射表（tab分割）
Generate sample-name to file-path mapping table (tab-separated)
"""

from .main import main
from .config import FofConfig
from .generator import FofGenerator

__version__ = "1.0.0"

__all__ = [
    'main',
    'FofConfig',
    'FofGenerator'
]
