"""候选基因RNA-seq转录验证模块|Candidate Gene RNA-seq Transcriptional Validation Module"""

from .config import GeneRnaseqCheckConfig
from .main import GeneRnaseqCheckPipeline

__version__ = "1.0.0"

__all__ = [
    'GeneRnaseqCheckConfig',
    'GeneRnaseqCheckPipeline',
]
