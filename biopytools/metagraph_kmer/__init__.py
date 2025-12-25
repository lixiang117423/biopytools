"""
🧬 MetaGraph K-mer库构建与查询分析工具包 | MetaGraph K-mer Library Construction and Query Analysis Toolkit

功能 | Features:
    - 使用MetaGraph构建k-mer图索引 | Build k-mer graph index using MetaGraph
    - 从FASTA文件保存坐标信息 | Save coordinate information from FASTA files
    - 使用KMC统计k-mer丰度 | Count k-mer abundance using KMC
    - 自动识别文件格式 | Automatic file format detection
    - Canonical k-mer处理 | Canonical k-mer processing
    
作者 | Author: Claude
版本 | Version: v2.0.0 - MetaGraph版本
日期 | Date: 2025-10-13

使用示例 | Usage Examples:
    from biopytools.metagraph_kmer import KmerAnalyzer, KmerConfig
    
    # 创建分析器 | Create analyzer
    analyzer = KmerAnalyzer(
        reference="reference.fasta",
        query="query.fastq",
        kmer_length=31,
        threads=88,
        output_dir="kmer_results"
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "2.0.0"
__author__ = "Claude"
__date__ = "2025-10-13"

from .main import KmerAnalyzer
from .config import KmerConfig

__all__ = ['KmerAnalyzer', 'KmerConfig']
