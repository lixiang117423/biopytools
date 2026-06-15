"""
LAI模块|LTR Assembly Index Module

用于评估基因组组装质量的长末端重复逆转座子组装指数计算工具
Long Terminal Repeat (LTR) Assembly Index calculator for genome assembly quality assessment

作者|Author: Xiang LI
版本|Version: 0.1.0

使用示例|Usage Examples:
    from biopytools.lai import LAICalculator, LAIConfig, LAILogger

    # 创建计算器|Create calculator
    calculator = LAICalculator(
        genome="genome.fa",
        output_dir="lai_output",
        threads=12
    )

    # 运行分析|Run analysis
    calculator.run()
"""

from .main import LAICalculator
from .config import (
    LAIConfig,
    LTRHarvestConfig,
    LTRFinderConfig,
    LTRRetrieverConfig,
    LAICalculateConfig
)
from .utils import LAILogger

__version__ = "0.1.0"
__author__ = "Xiang LI"

__all__ = [
    # 主类|Main class
    'LAICalculator',

    # 配置类|Configuration classes
    'LAIConfig',
    'LTRHarvestConfig',
    'LTRFinderConfig',
    'LTRRetrieverConfig',
    'LAICalculateConfig',

    # 日志器|Logger
    'LAILogger'
]
