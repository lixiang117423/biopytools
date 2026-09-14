"""
BLAST分析工具函数模块|BLAST Analysis Utility Functions Module
仅保留conda环境包装函数;样品映射/命令执行等逻辑由main.py的BLASTAnalyzer负责
|Only conda wrapping kept here; sample mapping / command execution live in BLASTAnalyzer (main.py)
"""

# conda包装统一走公共层(§13): 同源conda绝对路径 + run -p <环境前缀>, 严禁裸调conda
from ..common.conda_runner import build_conda_command
