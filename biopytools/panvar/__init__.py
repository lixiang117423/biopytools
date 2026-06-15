"""泛基因组变异分析模块|Pan-genome Variant Analysis Module

封装vg deconstruct + VCF变异分类统计 + 泛基因组增长曲线(gamma)完整流程
Wraps vg deconstruct + VCF variant classification + growth curve (gamma) pipeline

自动识别输入: .gbz执行deconstruct+统计, .vcf直接统计
Auto-detect input: .gbz runs deconstruct+summary, .vcf runs summary only

使用示例|Usage Examples:
    biopytools panvar -i cactus_output.gbz -P T2T -o output/
    biopytools panvar -i cactus_SV.vcf -o output/
"""

from .config import PanvarConfig
from .runner import PanvarRunner

__version__ = "1.0.0"
__author__ = "Xiang LI"

__all__ = ['PanvarConfig', 'PanvarRunner']
