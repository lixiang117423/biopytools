"""
Kmer矩阵比较工具包|Kmer Matrix Comparison Toolkit
功能: 比较两个kmer矩阵文件，找出特有kmer并计算窗口统计|
Features: Compare two kmer matrix files, find unique kmers and calculate window statistics
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-02

使用示例|Usage Examples:
    from biopytools.kmertools.compare import KmerMatrixComparator, KmerCompareConfig

    # 创建比较器|Create comparator
    comparator = KmerMatrixComparator(
        file1="matrix1.txt",
        file2="matrix2.txt",
        output_prefix="comparison",
        window_size=100000
    )

    # 运行比较|Run comparison
    comparator.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import KmerMatrixComparator
from .config import KmerCompareConfig

__all__ = ['KmerMatrixComparator', 'KmerCompareConfig']
