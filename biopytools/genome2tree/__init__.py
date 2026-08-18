"""genome2tree 模块|genome2tree module:基因组目录→物种树(waster)|genome dir to species tree"""

from .calculator import Genome2TreeCalculator
from .config import Genome2TreeConfig

__version__ = "1.0.0"
__all__ = ["Genome2TreeConfig", "Genome2TreeCalculator", "__version__"]
