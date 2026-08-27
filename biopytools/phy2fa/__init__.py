"""phy2fa 模块|phy2fa module:Phylip→FASTA 转换|Phylip to FASTA conversion"""

from .calculator import Phy2FaCalculator
from .config import Phy2FaConfig

__version__ = "1.0.0"
__all__ = ["Phy2FaConfig", "Phy2FaCalculator", "__version__"]
