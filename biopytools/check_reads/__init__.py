"""check_reads 模块|check_reads module:fastq 完整性检查|FASTQ integrity check"""

from .calculator import CheckReadsCalculator
from .config import CheckReadsConfig

__version__ = "1.0.0"
__all__ = ["CheckReadsConfig", "CheckReadsCalculator", "__version__"]
