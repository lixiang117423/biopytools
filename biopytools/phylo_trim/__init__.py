"""phylo_trim 整合工具包|phylo_trim Integration Toolkit
功能: 整合 mafft-fasttree + trimal,输出 trimal 前后两棵系统发育树
|Features: Integrate mafft-fasttree + trimal, output before/after-trimal trees
"""

__version__ = "1.0.0"

from .main import PhyloTrimRunner
from .config import PhyloTrimConfig

__all__ = [
    'PhyloTrimRunner',
    'PhyloTrimConfig',
]
