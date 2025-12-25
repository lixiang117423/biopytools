"""
Fastq到VCF (Parabricks) 模块 | Fastq to VCF (Parabricks) Module

基于Parabricks的重测序全基因组变异检测全流程分析模块
Whole genome resequencing variant detection pipeline based on Parabricks
"""

from .main import Fastq2VcfParabricksProcessor

__all__ = ['Fastq2VcfParabricksProcessor']
__version__ = '1.0.0'