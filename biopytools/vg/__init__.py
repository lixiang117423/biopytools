"""
VG变异图分析模块|VG Variation Graph Analysis Module

本模块提供了VG (Variation Graph) 工具的Python封装，支持子命令方式调用
This module provides Python wrapper for VG (Variation Graph) toolkit with subcommand support

主要功能|Main features:
    - construct: 从VCF和参考基因组构建变异图|Build variation graph from VCF and reference
    - giraffe: 快速序列比对|Fast read alignment
    - deconstruct: 从变异图生成VCF|Export VCF from variation graph
    - index: 创建图索引|Create graph indexes

作者|Author: Xiang LI
版本|Version: 1.0.0

使用示例|Usage Examples:
    # 构建变异图|Build variation graph
    biopytools vg construct -r ref.fa -v variants.vcf.gz -o graph.vg

    # 序列比对|Align reads
    biopytools vg giraffe -g graph.xg -f reads.fq -o alignments.gam

    # 导出VCF|Export VCF
    biopytools vg deconstruct -i graph.vg -r ref_path -o output.vcf
"""

from .main import main

__version__ = "1.0.0"
__all__ = ['main']
