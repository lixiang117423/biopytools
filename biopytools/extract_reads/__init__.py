"""
基于contig-reads对应关系从fastq提取reads|Extract reads from fastq by contig-reads mapping
功能: 根据contig与reads的对应关系表，从fastq文件中提取指定的reads并压缩输出|
Features: Extract specified reads from fastq file based on contig-reads mapping and compress output
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-01-10

使用示例|Usage Examples:
    from biopytools.extract_reads import ReadsExtractor, ExtractReadsConfig

    # 创建提取器|Create extractor
    extractor = ReadsExtractor(
        mapping_file="contig_reads.tsv",
        fastq_file="input.fq.gz",
        output_prefix="extracted"
    )

    # 运行提取|Run extraction
    extractor.extract()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import ReadsExtractor
from .config import ExtractReadsConfig
from .extractor import FastqReadsExtractor

__all__ = ['ReadsExtractor', 'ExtractReadsConfig', 'FastqReadsExtractor']
