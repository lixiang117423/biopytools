"""
Primer3引物设计工具包|Primer3 Primer Design Toolkit
功能: 批量PCR引物设计的完整流程|Features: Complete pipeline for batch PCR primer design
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-03-11

使用示例|Usage Examples:
    from biopytools.primer3 import Primer3Evaluator, Primer3Config

    # 创建评估器|Create evaluator
    evaluator = Primer3Evaluator(
        input_fasta="sequences.fasta",
        output_dir="./primer3_output"
    )

    # 运行引物设计|Run primer design
    evaluator.run_design()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .primer3_evaluator import Primer3Evaluator
from .config import Primer3Config
from .parser import FastaParser, Primer3InputGenerator, Primer3OutputParser, ResultsFormatter
from .utils import Primer3Logger, CommandRunner, get_conda_env, build_conda_command

__all__ = [
    'Primer3Evaluator',
    'Primer3Config',
    'FastaParser',
    'Primer3InputGenerator',
    'Primer3OutputParser',
    'ResultsFormatter',
    'Primer3Logger',
    'CommandRunner',
    'get_conda_env',
    'build_conda_command'
]
