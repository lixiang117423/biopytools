"""
覆盖度过滤工具包|Coverage Filter Toolkit
功能: 基于BAM覆盖度对序列进行质量分级和过滤|
Features: Sequence quality classification and filtering based on BAM coverage
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-01-26

使用示例|Usage Examples:
    from biopytools.coverage_filter import CoverageFilter, CoverageFilterConfig

    # 创建过滤器|Create filter
    filter_tool = CoverageFilter(
        bam_file="sample.bam",
        fasta_file="genome.fa",
        output_prefix="filtered"
    )

    # 运行过滤|Run filtering
    filter_tool.run_filter()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import CoverageFilter
from .config import CoverageFilterConfig

__all__ = ['CoverageFilter', 'CoverageFilterConfig']
