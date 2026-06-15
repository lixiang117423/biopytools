"""
WGDI工具包|WGDI Toolkit
功能: 比较基因组学分析工具封装|
Features: Comparative genomics analysis tools wrapper
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-03-09

使用示例|Usage Examples:
    from biopytools.wgdi import WGDIProcessor, DotPlotConfig

    # 创建DotPlot配置|Create DotPlot configuration
    config = DotPlotConfig(
        blast_file="blast.txt",
        gff1_file="species1.gff",
        gff2_file="species2.gff",
        lens1_file="species1.lens",
        lens2_file="species2.lens"
    )

    # 运行分析|Run analysis
    processor = WGDIProcessor(config)
    processor.run_dotplot(config)
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .processor import WGDIProcessor
from .config import (
    WGDIConfig,
    DotPlotConfig,
    CollinearityConfig,
    CalKsConfig
)
from .utils import WGDILogger, CommandRunner, WGDIConfGenerator

__all__ = [
    'WGDIProcessor',
    'WGDIConfig',
    'DotPlotConfig',
    'CollinearityConfig',
    'CalKsConfig',
    'WGDILogger',
    'CommandRunner',
    'WGDIConfGenerator'
]
