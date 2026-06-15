"""
染色体重命名工具包|Chromosome Rename Toolkit
功能: 基于minimap2比对结果进行染色体命名转换|
Features: Chromosome name conversion based on minimap2 alignment results
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-24

使用示例|Usage Examples:
    from biopytools.chr_rename import ChrRenameAnnotator, ChrRenameConfig

    # 创建重命名器|Create renamer
    renamer = ChrRenameAnnotator(
        ref_fasta="reference.fa",
        query_fasta="query.fa",
        output_dir="./output"
    )

    # 运行重命名|Run renaming
    renamer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import ChrRenameAnnotator
from .config import ChrRenameConfig

__all__ = ['ChrRenameAnnotator', 'ChrRenameConfig']
