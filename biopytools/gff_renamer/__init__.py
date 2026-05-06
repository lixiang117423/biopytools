"""
GFF文件重命名模块|GFF File Renamer Module
"""

from .main import main as gff_renamer_main
from .config import GFFRenamerConfig
from .renamer import GFFRenamer

__version__ = "1.0.0"

__all__ = [
    'gff_renamer_main',
    'GFFRenamerConfig',
    'GFFRenamer'
]
