"""
BRAKER3基因组注释工具包|BRAKER3 Genome Annotation Toolkit
功能: 基于BRAKER3的真核生物基因组基因结构注释完整流程|
Features: Complete pipeline for eukaryotic genome gene structure annotation based on BRAKER3
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-02

使用示例|Usage Examples:
    from biopytools.braker import BrakerPipeline, BrakerConfig

    # 创建注释流程|Create annotation pipeline
    pipeline = BrakerPipeline(
        genome="genome.fasta",
        species="my_oomycete",
        prot_seq="proteins.fasta",
        isoseq="isoseq.fasta",
        rnaseq_dirs=["/path/to/rnaseq1", "/path/to/rnaseq2"],
        threads=12
    )

    # 运行完整流程|Run complete pipeline
    pipeline.run_pipeline()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .pipeline import BrakerPipeline
from .config import BrakerConfig

__all__ = ['BrakerPipeline', 'BrakerConfig']
