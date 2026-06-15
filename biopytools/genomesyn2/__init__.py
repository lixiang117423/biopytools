"""
GenomeSyn2比较基因组学可视化工具包|GenomeSyn2 Comparative Genomics Visualization Toolkit
功能: 基因组共线性可视化、结构变异分析、祖先血统解析|
Features: Genome synteny visualization, structural variation analysis, ancestry deconvolution
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-07

使用示例|Usage Examples:
    from biopytools.genomesyn2 import GenomeSyn2Runner, GenomeSyn2Config

    # 创建运行器|Create runner
    runner = GenomeSyn2Runner(
        align='mummer',
        genome='./genome_dir/',
        outdir='./output/',
        threads=12
    )

    # 运行分析|Run analysis
    runner.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import GenomeSyn2Runner
from .config import GenomeSyn2Config
from .calculator import GenomeSyn2Calculator

__all__ = ['GenomeSyn2Runner', 'GenomeSyn2Config', 'GenomeSyn2Calculator']
