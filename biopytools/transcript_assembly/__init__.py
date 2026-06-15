"""
转录本从头组装工具包|Transcript De Novo Assembly Toolkit
功能: 基于HISAT2比对和StringTie组装的转录本从头组装流程|
Features: De novo transcript assembly pipeline based on HISAT2 alignment and StringTie assembly
作者|Author: Xiang LI
版本|Version: v1.0.0
日期|Date: 2026-04-25

使用示例|Usage Examples:
    from biopytools.transcript_assembly import TranscriptAssembler, TranscriptAssemblyConfig

    # 创建组装器|Create assembler
    assembler = TranscriptAssembler(
        genome_file="genome.fasta",
        input_dir="./clean_data",
        output_dir="./transcript_assembly_output"
    )

    # 运行组装|Run assembly
    assembler.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import TranscriptAssembler
from .config import TranscriptAssemblyConfig

__all__ = ['TranscriptAssembler', 'TranscriptAssemblyConfig']
