"""
Fastq到VCF (GTX) 模块 | Fastq to VCF (GTX) Module

基于GTX的重测序全基因组变异检测全流程分析模块
Whole genome resequencing variant detection pipeline based on GTX
"""

from .main import Fastq2VcfGTXProcessor

__all__ = ['Fastq2VcfGTXProcessor']
__version__ = '1.0.0'