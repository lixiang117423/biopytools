"""
HMMsearch结果处理工具包|HMMsearch Result Processing Toolkit
功能: 解析hmmsearch domtblout输出，生成易读表格并提取序列|
Features: Parse hmmsearch domtblout output, generate readable tables and extract sequences
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-28

使用示例|Usage Examples:
    from biopytools.hmmsearch import HMMsearchAnalyzer

    # 创建分析器|Create analyzer
    analyzer = HMMsearchAnalyzer(
        domtblout_file="results.domtblout",
        protein_fasta="proteins.fa",
        evalue_threshold=1e-10
    )

    # 运行分析|Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import HMMsearchAnalyzer
from .config import HMMsearchConfig
from .parser import DomtbloutParser, DomainHit
from .sequence_extractor import SequenceExtractor
from .hmmsearch_runner import HMMsearchRunner

__all__ = [
    'HMMsearchAnalyzer',
    'HMMsearchConfig',
    'DomtbloutParser',
    'DomainHit',
    'SequenceExtractor',
    'HMMsearchRunner'
]
