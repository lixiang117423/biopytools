"""
基因组组装质量评估模块|Genome Assembly Quality Control Module

提供基因组组装质量的综合评估，包括：
- BUSCO完整性评估
- LAI指数评估
- QV质量值计算
- Mapping评估
- 综合报告生成
"""

from .main import AssemblyQC, main

__version__ = "1.0.0"

__all__ = [
    "AssemblyQC",
    "main",
]
