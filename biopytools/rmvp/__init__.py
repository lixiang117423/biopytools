"""
rMVP GWAS分析模块|rMVP GWAS Analysis Module

支持GLM、MLM、FarmCPU三个模型的批量GWAS分析
Support batch GWAS analysis using GLM, MLM, and FarmCPU models

作者|Author: biopytools团队
版本|Version: 1.0.0
"""

from .config import RMVPConfig
from .analyzer import RMVPAnalyzer
from .result_parser import RMVPResultParser
from .cli import main

__version__ = "1.0.0"
__all__ = ["RMVPConfig", "RMVPAnalyzer", "RMVPResultParser", "main"]
