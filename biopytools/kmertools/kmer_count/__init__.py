"""
K-mer丰度分析工具包|K-mer Abundance Analysis Toolkit
功能: 使用Jellyfish统计k-mer丰度并生成存在/缺失矩阵|
Features: Count k-mer abundance using Jellyfish and generate presence/absence matrix
作者|Author: biopytools team
版本|Version: 1.0.0

使用示例|Usage Examples:
    from biopytools.kmertools.kmer_count import KmerCountAnalyzer, KmerCountConfig

    # 创建配置|Create configuration
    config = KmerCountConfig(
        input_dir="~/fastq_dir",
        pattern="*_1.fq.gz",
        kmer_lib="kmers.fasta",
        output_dir="results/"
    )

    # 创建分析器并运行|Create analyzer and run
    analyzer = KmerCountAnalyzer(config)
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "biopytools team"

from .main import KmerCountAnalyzer, main
from .config import KmerCountConfig
from .utils import KmerCountLogger

__all__ = [
    'KmerCountAnalyzer',
    'KmerCountConfig',
    'KmerCountLogger',
    'main'
]
