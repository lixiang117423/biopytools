"""
KMC k-mer分析工具包|KMC K-mer Analysis Toolkit
功能: 基于KMC的k-mer统计、丰度矩阵构建和查询|
Features: K-mer counting, abundance matrix construction and query based on KMC
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2025-01-22

使用示例|Usage Examples:
    from biopytools.kmc import KMCCounter, KMCConfig, KMCMatrixBuilder, KMCQuery

    # 配置|Configuration
    config = KMCConfig(
        input_files=["sample1.fq", "sample2.fq"],
        kmer_size=21,
        output_dir="./kmc_output"
    )

    # 统计k-mer|Count k-mers
    counter = KMCCounter(config)
    counter.run()

    # 构建丰度矩阵|Build abundance matrix
    matrix_builder = KMCMatrixBuilder(config)
    matrix_builder.build_matrix()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import KMCManager
from .config import KMCConfig
from .kmer_counter import KMCCounter
from .matrix_builder import KMCMatrixBuilder
from .kmer_query import KMCQuery

__all__ = [
    'KMCManager',
    'KMCConfig',
    'KMCCounter',
    'KMCMatrixBuilder',
    'KMCQuery'
]
