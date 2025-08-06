"""
序列提取工具包 | Sequence Extraction Toolkit
功能: 从VCF文件和基因组文件中提取特定区间的序列变异信息 | 
Features: Extract sequence variation information from VCF and genome files for specific regions
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-07-21

使用示例 | Usage Examples:
    from biopytools.sequence_toolkit import SequenceExtractor, SequenceConfig
    
    # 创建提取器 | Create extractor
    extractor = SequenceExtractor(
        vcf_file="variants.vcf",
        genome_file="genome.fa",
        chrom="chr1",
        start=1000,
        end=1050,
        output_dir="results"
    )
    
    # 运行提取 | Run extraction
    extractor.run_extraction()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import SequenceExtractor
from .config import SequenceConfig

__all__ = ['SequenceExtractor', 'SequenceConfig']
