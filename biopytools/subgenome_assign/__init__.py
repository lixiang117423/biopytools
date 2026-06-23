"""
亚基因组归属工具包|Subgenome Assignment Toolkit

通过将目标多倍体基因组的每条染色体比对到各亲本参考，
按比对覆盖度判定每条染色体的亚基因组来源
|Assign each chromosome of a target polyploid genome to a subgenome
by alignment against parental references
"""

from .main import SubgenomeAssignRunner
from .config import SubgenomeAssignConfig

__version__ = "1.0.0"

__all__ = [
    'SubgenomeAssignRunner',
    'SubgenomeAssignConfig',
]
