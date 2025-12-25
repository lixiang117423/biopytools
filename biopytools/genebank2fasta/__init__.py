"""
GenBank序列提取工具包 🧬 | GenBank Sequence Extraction Toolkit
功能: 从GenBank文件批量提取CDS序列和蛋白质序列，支持多样品处理和统计分析 | 
Features: Batch extraction of CDS and protein sequences from GenBank files, supporting multi-sample processing and statistical analysis
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-09-26

使用示例 | Usage Examples:
    from biopytools.genebank2fasta import GenBankExtractor, ExtractorConfig
    
    # 创建提取器 | Create extractor
    extractor = GenBankExtractor(
        input_dir="/path/to/genbank/files",
        output_dir="/path/to/output",
        threads=88
    )
    
    # 运行提取 | Run extraction
    extractor.run_extraction()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import GenBankExtractor
from .config import ExtractorConfig

__all__ = ['GenBankExtractor', 'ExtractorConfig']
