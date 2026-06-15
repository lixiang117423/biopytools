"""
序列到基因组比对模块|Sequence to Genome Alignment Module

自动检测序列类型（DNA或蛋白质）并选择合适的比对工具：
- DNA序列：使用Minimap2进行比对
- 蛋白质序列：使用Miniprot进行比对
支持PAF格式解析、统计分析和GFF3/BED格式导出

Auto-detect sequence type (DNA or protein) and select appropriate alignment tool:
- DNA sequences: Align using Minimap2
- Protein sequences: Align using Miniprot
With PAF parsing, statistics analysis, and GFF3/BED export

作者|Author: Biopytools Development Team
版本|Version: 2.0.0
"""

__version__ = "2.0.0"
__author__ = "Biopytools Development Team"

from .main import Seq2GenomeAnalyzer, main
from .config import Seq2GenomeConfig
from .parser import PAFParser, PAFRecord
from .calculator import MiniprotAligner, Minimap2Aligner
from .statistics import StatisticsGenerator
from .gff3_exporter import GFF3Exporter, BEDExporter
from .utils import detect_sequence_type

__all__ = [
    "Seq2GenomeAnalyzer",
    "main",
    "Seq2GenomeConfig",
    "PAFParser",
    "PAFRecord",
    "MiniprotAligner",
    "Minimap2Aligner",
    "StatisticsGenerator",
    "GFF3Exporter",
    "BEDExporter",
    "detect_sequence_type",
]
