"""
HiCanu基因组组装工具包|HiCanu Genome Assembly Toolkit
功能: 使用Canu进行基因组组装，特别是PacBio HiFi reads
Features: Genome assembly using Canu, especially for PacBio HiFi reads

作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2025-12-30

使用示例|Usage Examples:
    from biopytools.hicanu import HiCanuPipeline, HiCanuConfig

    # 创建组装流程|Create assembly pipeline
    pipeline = HiCanuPipeline(
        reads_file="reads.fastq",
        genome_size="120m",
        prefix="sample1",
        output_dir="./hicanu_output"
    )

    # 运行组装|Run assembly
    pipeline.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import HiCanuPipeline
from .config import HiCanuConfig
from .calculator import HiCanuCalculator
from .utils import CanuLogger

__all__ = [
    'HiCanuPipeline',
    'HiCanuConfig',
    'HiCanuCalculator',
    'CanuLogger'
]
