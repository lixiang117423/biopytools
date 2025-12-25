"""
GFF格式转换工具包 | GFF Format Converter Toolkit
功能: 标准化GFF文件格式，重新编号基因ID并简化属性结构 | 
Features: Standardize GFF file format, renumber gene IDs and simplify attribute structure
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-09-04

使用示例 | Usage Examples:
    from biopytools.gff_converter import GFFConverter, GFFConfig
    
    # 创建转换器 | Create converter
    converter = GFFConverter(
        input_file="input.gff",
        output_file="output.gff",
        species_name="OV53",
        species_prefix="Ov"
    )
    
    # 运行转换 | Run conversion
    converter.run_conversion()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import GFFConverter
from .config import GFFConfig

__all__ = ['GFFConverter', 'GFFConfig']
