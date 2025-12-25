"""
HapHiC基因组scaffolding工具模块 | HapHiC Genome Scaffolding Tool Module

基于Hi-C数据的快速、参考基因组独立的等位基因感知scaffolding工具。
支持单倍型分相基因组组装、二倍体和多倍体基因组组装。
Fast, reference-independent, allele-aware scaffolding tool based on Hi-C data.
Supports haplotype-phased, haplotype-collapsed diploid and allopolyploid genome assemblies.
"""

__version__ = "1.0.0"
__author__ = "biopytools development team"

from .main import HapHiCProcessor
from .config import HapHiCConfig

__all__ = ['HapHiCProcessor', 'HapHiCConfig']