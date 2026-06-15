"""
启动子提取器模块|Promoter Extractor Module

从GFF3注释文件和基因组FASTA文件中提取基因启动子序列|Extract gene promoter sequences from GFF3 annotation and genome FASTA files

功能特性|Features:
- 支持GFF3格式基因注释|Support GFF3 format gene annotation
- 自动处理基因方向（+/-链）|Automatically handle gene strand (+/-)
- 边界自动截断处理|Automatic boundary truncation
- 输出FASTA序列和BED位置文件|Output FASTA sequences and BED position files
- 支持指定基因列表提取|Support extracting specified gene list
"""

from .main import PromoterRunner, main
from .config import PromoterExtractorConfig
from .calculator import PromoterExtractor, GFF3Parser, GenomeSequenceLoader, GeneInfo

__version__ = "1.0.0"

__all__ = [
    'PromoterRunner',
    'main',
    'PromoterExtractorConfig',
    'PromoterExtractor',
    'GFF3Parser',
    'GenomeSequenceLoader',
    'GeneInfo'
]
