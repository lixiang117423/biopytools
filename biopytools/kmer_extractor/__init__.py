"""
🧬 K-mer提取工具包 | K-mer Extraction Toolkit
功能: 从FASTA/FASTQ文件中提取k-mer序列，支持反向互补和多文件批处理 | 
Features: Extract k-mer sequences from FASTA/FASTQ files with reverse complement and batch processing support
作者 | Author: Claude  
版本 | Version: v1.1 - 重构版：FASTA使用unikmer，FASTQ使用jellyfish | Refactored version: unikmer for FASTA, jellyfish for FASTQ
日期 | Date: 2025-08-11

使用示例 | Usage Examples:
    from biopytools.kmer_extractor import KmerExtractor, KmerConfig
    
    # 创建提取器 | Create extractor
    extractor = KmerExtractor(
        input_files=["sample1.fastq", "sample2.fastq"],
        output_dir="kmer_results",
        kmer_length=51,
        threads=88
    )
    
    # 运行提取 | Run extraction
    extractor.run_extraction()
"""

__version__ = "1.1.0"
__author__ = "Claude"

from .main import KmerExtractor
from .config import KmerConfig

__all__ = ['KmerExtractor', 'KmerConfig']
