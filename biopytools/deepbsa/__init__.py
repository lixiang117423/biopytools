"""
DeepBSA批量分析工具模块|DeepBSA Batch Analysis Tool Module

提供DeepBSA多种BSA分析方法的批量运行功能
Provides batch running functionality for multiple DeepBSA BSA analysis methods

子命令|Sub-commands:
  - batch: 生成批量处理命令|Generate batch processing commands
  - run: 运行DeepBSA分析方法|Run DeepBSA analysis methods
  - merge: 合并分析结果|Merge analysis results
"""

from .main import main
from .cli import main as cli_main

__version__ = "2.0.0"

__all__ = ['main', 'cli_main']
