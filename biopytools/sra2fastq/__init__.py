"""
SRA转FASTQ转换工具包|SRA to FASTQ Conversion Toolkit
功能: 使用parallel-fastq-dump高速批量转换SRA文件，支持真正的多线程加速|Features: High-speed batch conversion using parallel-fastq-dump with true multi-threading
作者|Author: Claude
版本|Version: 2.0.0
日期|Date: 2025-10-11

使用示例|Usage Examples:
    from biopytools.sra_converter import SRAConverter, ConvertConfig
    
    # 创建转换器|Create converter
    converter = SRAConverter(
        input_path="sra_files/",
        output_dir="fastq_output",
        compress=True,
        threads=12
    )
    
    # 运行转换|Run conversion
    converter.run_conversion()
"""

__version__ = "2.0.0"
__author__ = "Claude"

from .main import SRAConverter
from .config import ConvertConfig

__all__ = ['SRAConverter', 'ConvertConfig']
