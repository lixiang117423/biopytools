"""
VCF抽样工具包|VCF Sampling Toolkit
功能: 从VCF文件中按比例随机抽取SNP位点 |
Features: Randomly sample SNP sites from VCF files by proportion
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2025-01-03

使用示例|Usage Examples:
    from biopytools.vcf_sampler import VCFSampler, VCFSamplerConfig

    # 创建抽样器|Create sampler
    sampler = VCFSampler(
        input_vcf="input.vcf.gz",
        output_vcf="output.vcf.gz",
        sample_rate=0.1
    )

    # 运行抽样|Run sampling
    sampler.run_sampling()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import VCFSampler
from .config import VCFSamplerConfig

__all__ = ['VCFSampler', 'VCFSamplerConfig']
