"""
VCF转换工具包|VCF Converter Toolkit
功能: 将VCF格式的SNP数据转换为多种系统发生学矩阵格式|
Features: Convert VCF format SNP data to various phylogenetic matrix formats
作者|Author: Based on Edgardo M. Ortiz's script, modularized version
版本|Version: v2.9.1 - 模块化版本|Modular version
日期|Date: 2025-08-25

支持格式|Supported formats:
- PHYLIP: 系统发生学分析标准格式
- FASTA: 序列分析通用格式  
- NEXUS: 包含元数据的系统发生学格式
- Binary NEXUS: 适用于SNAPP的二进制格式

使用示例|Usage Examples:
    from vcf_converter import VCFConverter
    
    # 创建转换器|Create converter
    converter = VCFConverter(
        input_file="variants.vcf.gz",
        output_dir="converted_results",
        threads=12
    )
    
    # 运行转换|Run conversion
    converter.run_conversion()
"""

__version__ = "2.9.1"
__author__ = "Based on Edgardo M. Ortiz's script"

from .main import VCFConverter
from .config import ConverterConfig

__all__ = ['VCFConverter', 'ConverterConfig']
