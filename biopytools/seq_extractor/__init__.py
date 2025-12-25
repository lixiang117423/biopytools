"""
序列提取工具包 🧬 | Sequence Extraction Toolkit
功能: 从DNA或蛋白质序列文件中提取指定区域序列的完整工具 | 
Features: Complete tool for extracting specified regions from DNA or protein sequence files
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-09-27

使用示例 | Usage Examples:
    from biopytools.seq_extractor import SequenceExtractor, ExtractorConfig
    
    # 创建提取器 | Create extractor
    extractor = SequenceExtractor(
        sequence_file="genome.fasta",
        regions_file="regions.bed",
        output_file="extracted_sequences.fasta",
        sequence_type="dna"
    )
    
    # 运行提取 | Run extraction
    extractor.run_extraction()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import SequenceExtractor
from .config import ExtractorConfig

__all__ = ['SequenceExtractor', 'ExtractorConfig']
