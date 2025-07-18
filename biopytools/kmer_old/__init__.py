"""
高性能k-mer数据库分析工具包 | High-Performance K-mer Database Analysis Toolkit
基于kmtricks + RocksDB的大规模k-mer分析系统 | Large-scale k-mer analysis system based on kmtricks + RocksDB
适用于大规模数据集（数千个样本）| Suitable for large-scale datasets (thousands of samples)

使用示例 | Usage Examples:
    from kmer import KmerAnalyzer, KmerConfig
    
    # 创建分析器 | Create analyzer
    analyzer = KmerAnalyzer(
        gene_fasta="genes.fasta",
        fastq_dir="/path/to/fastq",
        output_dir="./results",
        kmer_size=51
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "wheatomics"

from .main import KmerAnalyzer
from .config import KmerConfig

__all__ = ['KmerAnalyzer', 'KmerConfig']
