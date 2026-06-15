"""
FASTQ文件过滤工具包|FASTQ File Filtering Toolkit
功能: 基于GC含量和序列长度过滤FASTQ文件|Features: Filter FASTQ files by GC content and sequence length
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-03

使用示例|Usage Examples:
    from biopytools.fastq_gc_filter import FastqGcFilter

    # 创建过滤器|Create filter
    filter_tool = FastqGcFilter(
        input_file="reads.fastq",
        output_file="filtered.fastq",
        min_gc=25.0,
        max_gc=70.0,
        min_length=50
    )

    # 运行过滤|Run filtering
    filter_tool.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import FastqGcFilter
from .config import FastqGcFilterConfig

__all__ = ['FastqGcFilter', 'FastqGcFilterConfig']
