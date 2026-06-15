"""
Pixy群体遗传学统计模块|Pixy Population Genetics Statistics Module

提供pi、dxy、fst等群体遗传学统计量的计算功能
Provides calculation of population genetics statistics such as pi, dxy, fst

作者|Author: Xiang LI
版本|Version: 1.0.0

使用示例|Usage Examples:
    from biopytools.pixy import main, PixyConfig, PixyLogger

    # 通过命令行使用|Use via command line
    # biopytools pixy -i variants.vcf.gz -p populations.txt -o pixy_output
"""

from .main import main
from .config import PixyConfig
from .utils import PixyLogger

__version__ = "1.0.0"
__author__ = "Xiang LI"

__all__ = [
    'main',
    'PixyConfig',
    'PixyLogger'
]
