"""
VCF转系统发育树工具包|VCF to Phylogenetic Tree Toolkit

功能|Features:
  从VCF SNP数据直接构建系统发育树，支持FastTree和IQ-TREE后端
  Direct phylogenetic tree construction from VCF SNP data, supporting FastTree and IQ-TREE backends

作者|Author: Claude
版本|Version: v1.0.0
日期|Date: 2026-08-07

使用示例|Usage Examples:
    from biopytools.vcf2tree import Vcf2TreeRunner, Vcf2TreeConfig

    # 默认IQ-TREE建树|Default IQ-TREE
    runner = Vcf2TreeRunner(
        input_file="variants.vcf",
        output_dir="tree_output",
        threads=12
    )
    runner.run_pipeline()

    # FastTree建树|FastTree
    runner = Vcf2TreeRunner(
        input_file="variants.vcf",
        output_dir="tree_output",
        method="fasttree",
        threads=12
    )
    runner.run_pipeline()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import Vcf2TreeRunner
from .config import Vcf2TreeConfig
from .utils import Vcf2TreeLogger

__all__ = [
    'Vcf2TreeRunner',
    'Vcf2TreeConfig',
    'Vcf2TreeLogger'
]
