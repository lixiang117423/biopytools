"""
Wgsim Module
Wgsim基因组测序数据模拟模块

功能: 使用wgsim从基因组参考序列模拟Illumina双端测序reads

使用示例|Usage Examples:
    from biopytools.wgsim import WgsimRunner

    runner = WgsimRunner(
        input_dir="01.data/genome",
        output_dir="01.data/raw"
    )
    runner.run()

作者|Author: Xiang Li
版本|Version: 1.0.0
日期|Date: 2026-06-02
"""

__version__ = "1.0.0"
__author__ = "Xiang Li"

from .main import WgsimRunner
from .config import WgsimConfig

__all__ = ['WgsimRunner', 'WgsimConfig']
