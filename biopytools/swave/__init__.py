"""
Swave结构变异检测工具包|Swave Structural Variant Detection Toolkit
功能: 从泛基因组图检测结构变异和复杂SV|
Features: Detect structural variants and complex SVs from pangenome graphs
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-03-16

使用示例|Usage Examples:
    from biopytools.swave import SwaveRunner, SwaveConfig

    # 创建SV检测器|Create SV caller
    runner = SwaveRunner(
        assemblies_tsv="assemblies.tsv",
        ref_fasta="reference.fa",
        gfa_file="pangenome.gfa",
        gfa_source="minigraph"
    )

    # 运行检测|Run detection
    runner.run_call()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import SwaveRunner
from .config import SwaveConfig
from .sv_caller import SwaveSVCaller, SwaveConverter
from .pav_extractor import PAVExtractor

__all__ = ['SwaveRunner', 'SwaveConfig', 'SwaveSVCaller', 'SwaveConverter', 'PAVExtractor']
