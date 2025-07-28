"""
VCF单体型提取工具包 | VCF Haplotype Extractor Toolkit
功能: 从VCF文件提取指定位点的单体型信息，支持位置文件和标准化输出格式 | 
Features: Extract haplotype information from VCF files for specified loci, supporting position files and standardized output format
作者 | Author: Claude  
版本 | Version: v2.0 - 模块化版本 | Modular version
日期 | Date: 2025-07-28

使用示例 | Usage Examples:
    from biopytools.haplotype_extractor import HaplotypeExtractor, HaplotypeConfig
    
    # 使用位置文件 | Using position file
    extractor = HaplotypeExtractor(
        vcf_file="variants.vcf.gz",
        position_file="positions.txt",
        output_file="haplotypes.txt"
    )
    
    # 单个位点 | Single position
    extractor = HaplotypeExtractor(
        vcf_file="variants.vcf.gz",
        chromosome="OV12",
        position=93635286,
        output_file="haplotypes.txt"
    )
    
    # 运行分析 | Run analysis
    extractor.run_analysis()

位置文件格式 | Position File Format:
    # 可以有表头（推荐）| Can have header (recommended)
    CHR	POS
    OV12	93635286
    OV12	93440535
    OV12	93234780
    
    # 也可以没有表头 | Can also without header
    OV12	93635286
    OV12	93440535
    OV12	93234780

输出格式 | Output Format:
    CHROM	POS	REF	ALT	Sample1	Sample2	Sample3
    OV12	93635286	C	A	CA	AA	AA
    OV12	93440535	C	A	CA	AA	CA
"""

__version__ = "2.0.0"
__author__ = "Claude"

from .main import HaplotypeExtractor
from .config import HaplotypeConfig

__all__ = ['HaplotypeExtractor', 'HaplotypeConfig']
