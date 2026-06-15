"""
K-mer工具模块|K-mer Tools Module

提供k-mer索引构建、查询和分析功能
Provides k-mer indexing, query and analysis functions

作者|Author: Xiang LI
版本|Version: 1.0.0

使用示例|Usage Examples:
    from biopytools.kmertools import main, BuildConfig, QueryConfig, KmerToolsLogger

    # 通过命令行使用|Use via command line
    # 构建 k-mer 库|Build k-mer database
    # kmertools build -i ./fastq_dir -o ./kmer_db

    # 查询 k-mer 库|Query k-mer database
    # kmertools query -d ./kmer_db/rocksdb -q query.fa -o ./results
"""

from .main import main
from .config import (
    KmerToolsConfig,
    BuildConfig,
    QueryConfig,
    SplitFastaConfig,
    GenFofConfig,
    ImportDBConfig,
    ExtractConfig
)
from .utils import KmerToolsLogger

__version__ = "1.0.0"
__author__ = "Xiang LI"

__all__ = [
    # 主函数|Main function
    'main',

    # 配置类|Configuration classes
    'KmerToolsConfig',
    'BuildConfig',
    'QueryConfig',
    'SplitFastaConfig',
    'GenFofConfig',
    'ImportDBConfig',
    'ExtractConfig',

    # 日志器|Logger
    'KmerToolsLogger'
]
