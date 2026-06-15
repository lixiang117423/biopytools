"""
Minigraph泛基因组图构建和分析工具包|Minigraph Pangenome Graph Construction and Analysis Toolkit
功能: 构建泛基因组图、序列映射、SV调用和bubble提取|
Features: Build pangenome graphs, sequence mapping, SV calling and bubble extraction
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-03-16

使用示例|Usage Examples:
    from biopytools.minigraph import MinigraphRunner

    # 构建泛基因组图|Build pangenome graph
    runner = MinigraphRunner(
        command='build',
        ref_fasta='reference.fa',
        sample_fastas=['sample1.fa', 'sample2.fa'],
        output_gfa='pangenome.gfa'
    )
    runner.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import MinigraphRunner
from .config import (MinigraphBuildConfig, MinigraphCallConfig,
                     MinigraphBubbleConfig, MinigraphMapConfig)
from .graph_builder import (MinigraphGraphBuilder, MinigraphSVCaller,
                             MinigraphBubbleExtractor, MinigraphMapper)

__all__ = [
    'MinigraphRunner',
    'MinigraphBuildConfig',
    'MinigraphCallConfig',
    'MinigraphBubbleConfig',
    'MinigraphMapConfig',
    'MinigraphGraphBuilder',
    'MinigraphSVCaller',
    'MinigraphBubbleExtractor',
    'MinigraphMapper',
]
