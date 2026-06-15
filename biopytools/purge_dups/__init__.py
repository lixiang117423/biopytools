"""
Purge_Dups基因组去冗余工具包|Purge_Dups Genome Deduplication Toolkit
功能: 基于测序深度去除基因组组装中的单倍型和重叠序列|
Features: Remove haplotigs and overlaps in genome assembly based on read depth
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-01-21

使用示例|Usage Examples:
    from biopytools.purge_dups import PurgeDupsRunner, PurgeDupsConfig

    # 创建去冗余器|Create deduplicator
    runner = PurgeDupsRunner(
        input="assembly.fa",
        reads="pacbio_reads.fq",
        output_dir="purge_dups_output"
    )

    # 运行完整流程|Run complete pipeline
    runner.run_full_pipeline()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import PurgeDupsRunner
from .config import PurgeDupsConfig

__all__ = ['PurgeDupsRunner', 'PurgeDupsConfig']
