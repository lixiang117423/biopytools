"""
BAM to FASTQ转换模块|BAM to FASTQ Conversion Module

使用bam2fastq将BAM文件批量转换为FASTQ格式
Batch convert BAM files to FASTQ format using bam2fastq
"""

from .main import main
from .config import BAM2FASTQConfig
from .converter import BAMConverter

__version__ = "2.0.0"

__all__ = [
    'main',
    'BAM2FASTQConfig',
    'BAMConverter'
]
