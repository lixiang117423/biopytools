"""
蛋白质到基因组比对模块|Protein to Genome Alignment Module

使用Miniprot进行蛋白质与基因组的比对，支持PAF格式解析、统计分析和GFF3/BED格式导出
Align proteins to genome using Miniprot, with PAF parsing, statistics analysis, and GFF3/BED export

作者|Author: Biopytools Development Team
版本|Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "Biopytools Development Team"

from .main import Pep2GenomeAnalyzer, main
from .config import Pep2GenomeConfig
from .parser import PAFParser, PAFRecord
from .calculator import MiniprotAligner
from .statistics import StatisticsGenerator
from .gff3_exporter import GFF3Exporter, BEDExporter

__all__ = [
    "Pep2GenomeAnalyzer",
    "main",
    "Pep2GenomeConfig",
    "PAFParser",
    "PAFRecord",
    "MiniprotAligner",
    "StatisticsGenerator",
    "GFF3Exporter",
    "BEDExporter",
]
