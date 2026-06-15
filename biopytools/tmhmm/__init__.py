"""
TMHMM Module
TMHMM跨膜螺旋预测模块

功能: 使用TMHMM预测蛋白质序列的跨膜螺旋结构

使用示例|Usage Examples:
    from biopytools.tmhmm import TmhmmRunner

    runner = TmhmmRunner(
        input_file="proteins.fa",
        output_dir="output"
    )
    runner.run()

作者|Author: Xiang Li
版本|Version: 1.0.0
日期|Date: 2026-06-04
"""

__version__ = "1.0.0"
__author__ = "Xiang Li"

from .main import TmhmmRunner
from .config import TmhmmConfig

__all__ = ['TmhmmRunner', 'TmhmmConfig']
