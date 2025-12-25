"""
Dsuite Module
Dsuite D统计量分析模块

功能: 使用Dsuite进行基因渗入分析，计算D统计量(ABBA-BABA test)
      检测群体间的基因渗入事件

使用示例 | Usage Examples:
    from biopytools.dsuite import DsuiteAnalyzer

    analyzer = DsuiteAnalyzer(
        vcf_file="variants.vcf.gz",
        sets_file="sets.txt",
        output_dir="output"
    )
    analyzer.run()

作者 | Author: Xiang Li
版本 | Version: 1.0.0
日期 | Date: 2025-12-25
"""

__version__ = "1.0.0"
__author__ = "Xiang Li"

from .main import DsuiteAnalyzer
from .config import DsuiteConfig

__all__ = ['DsuiteAnalyzer', 'DsuiteConfig']
