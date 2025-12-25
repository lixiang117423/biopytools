"""
Rename Chromosomes Module
染色体重命名模块

功能: 将FASTA文件的序列重命名，前N条命名为Chr01, Chr02...
      剩余序列命名为HiC_scaffold_01, HiC_scaffold_02...

使用示例 | Usage Examples:
    from biopytools.rename_chromosomes import ChromosomeRenamer

    renamer = ChromosomeRenamer(
        input_file="genome.fa",
        output_file="renamed.fa",
        chromosome_number=20
    )
    renamer.rename()

作者 | Author: Xiang Li
版本 | Version: 1.0.0
日期 | Date: 2025-12-25
"""

__version__ = "1.0.0"
__author__ = "Xiang Li"

from .main import ChromosomeRenamer
from .config import ChromosomeRenameConfig

__all__ = ['ChromosomeRenamer', 'ChromosomeRenameConfig']
