"""
VCF按染色体合并工具包|VCF Merge by Chromosome Toolkit
功能: 自动按染色体合并分割的VCF文件|Features: Automatically merge split VCF files by chromosome
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2025-12-30

使用示例|Usage Examples:
    from biopytools.vcf_merger import VCFMerger, VCFMergerConfig

    # 创建合并器|Create merger
    merger = VCFMerger(
        input_dir="/path/to/vcf_files",
        output_dir="/path/to/output",
        pattern="*.joint.vcf.gz",
        threads=12
    )

    # 运行合并|Run merge
    merger.run()

    # 或者只获取分组信息|Or only get grouping info
    # merger.group_vcf_files()
    # chr_groups = merger.get_chromosome_groups()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import VCFMerger
from .config import VCFMergerConfig

__all__ = ['VCFMerger', 'VCFMergerConfig']
