"""序列提取工具(seqkit封装)|Sequence extraction tool (seqkit wrapper)"""

from .config import SeqExtractConfig
from .utils import SeqExtractRunner

__version__ = "1.0.0"

__all__ = [
    "SeqExtractConfig",
    "SeqExtractRunner",
]
