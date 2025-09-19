"""
最长转录本提取工具包 | Longest mRNA Transcript Extraction Toolkit
功能: 从基因组和GFF3注释文件中提取最长的mRNA转录本序列 | 
Features: Extract longest mRNA transcript sequences from genome and GFF3 annotation files
作者 | Author: Xiang LI  
版本 | Version: v1.0 - 模块化重构版 | Modular refactored version
日期 | Date: 2025-07-16

使用示例 | Usage Examples:
    from biopytools.longest_mrna import LongestMRNAExtractor, LongestMRNAConfig
    
    # 创建提取器 | Create extractor
    extractor = LongestMRNAExtractor(
        genome_file="genome.fa",
        gff3_file="annotation.gff3",
        output_file="longest_transcripts.fa"
    )
    
    # 运行提取 | Run extraction
    extractor.run_extraction()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import LongestMRNAExtractor
from .config import LongestMRNAConfig

__all__ = ['LongestMRNAExtractor', 'LongestMRNAConfig']
