"""
AGP转表格工具包|AGP to Table Converter Toolkit
功能: 将AGP格式文件转换为易读的表格格式|Features: Convert AGP format files to readable table format
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-03-20

使用示例|Usage Examples:
    from biopytools.agp2table import AGPConverter, AGPConfig

    # 创建转换器|Create converter
    converter = AGPConverter(
        agp_file="assembly.agp",
        output_file="assembly_table.txt"
    )

    # 运行转换|Run conversion
    converter.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import AGPConverter
from .config import AGPConfig

__all__ = ['AGPConverter', 'AGPConfig']
