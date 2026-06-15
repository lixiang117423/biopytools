"""
FASTQ配对修复工具包|FASTQ Pair Fixing Toolkit
功能: 批量修复配对混乱的FASTQ文件|Features: Batch fix paired-end FASTQ files with pairing issues
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-04-04

使用示例|Usage Examples:
    from biopytools.pair_fastq import FastqPairFixer, FastqPairConfig

    # 创建配对修复器|Create pair fixer
    fixer = FastqPairFixer(
        input_dir="raw_data",
        output_dir="fixed_data",
        threads=16
    )

    # 运行配对修复|Run pair fixing
    fixer.process_all_samples()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import FastqPairFixer
from .config import FastqPairConfig

__all__ = ['FastqPairFixer', 'FastqPairConfig']
