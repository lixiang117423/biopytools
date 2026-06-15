"""
JCVI共线性分析工具集|JCVI Synteny Analysis Toolkit

提供MCscan共线性分析、等位基因鉴定、宏观/微观共线性可视化功能
Provides MCscan collinearity analysis, allelic gene identification,
and macro/micro synteny visualization

Usage:
    biopytools jcvi mcscan  -i input_dir -o output_dir
    biopytools jcvi allelic -i input_dir -o output_dir
    biopytools jcvi macro   -i pair_dir --pairs A,B
    biopytools jcvi micro   -i pair_dir --gene XXX
"""

from .config import JcviBaseConfig
from .utils import (
    JcviLogger,
    build_jcvi_command,
    build_conda_command,
    discover_samples,
    get_sample_name,
)

__version__ = "2.0.0"
__all__ = [
    "JcviBaseConfig", "JcviLogger",
    "build_jcvi_command", "build_conda_command",
    "discover_samples", "get_sample_name",
]
