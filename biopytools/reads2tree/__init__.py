"""reads2tree 模块|reads2tree module:fastq 目录→物种树(WASTER)|FASTQ dir to species tree"""

from .calculator import Reads2TreeCalculator
from .config import Reads2TreeConfig

__version__ = "1.0.0"
__all__ = ["Reads2TreeConfig", "Reads2TreeCalculator", "__version__"]
