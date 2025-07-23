"""
FASTQ质控工具包 | FASTQ Quality Control Toolkit
功能: 使用fastp进行FASTQ文件批量质控的完整流程 | 
Features: Complete pipeline for batch FASTQ quality control using fastp
作者 | Author: Xiang LI  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-07-15

使用示例 | Usage Examples:
    from fastp import FastpProcessor, FastpConfig
    
    # 创建处理器 | Create processor
    processor = FastpProcessor(
        input_dir="./raw_data",
        output_dir="./clean_data",
        fastp_path="/path/to/fastp",
        threads=12
    )
    
    # 运行批处理 | Run batch processing
    processor.run_batch_processing()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import FastpProcessor
from .config import FastpConfig

__all__ = ['FastpProcessor', 'FastpConfig']
