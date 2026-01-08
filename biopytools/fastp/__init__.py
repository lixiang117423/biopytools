"""
FASTQ质控工具包|FASTQ Quality Control Toolkit

功能: 使用fastp进行FASTQ文件批量质控的完整流程
Features: Complete pipeline for batch FASTQ quality control using fastp

作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2024-12-30

使用示例|Usage Examples:
    from fastp import FastpProcessor, FastpConfig

    # 创建处理器|Create processor
    processor = FastpProcessor(
        input_dir="./raw_data",
        output_dir="./clean_data",
        fastp_path="/path/to/fastp",
        threads=12
    )

    # 运行批处理|Run batch processing
    processor.run_batch_processing()

命令行使用|Command Line Usage:
    python -m fastp.main -i raw_data/ -o clean_data/ -t 12
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import FastpProcessor
from .config import FastpConfig

__all__ = ['FastpProcessor', 'FastpConfig']
