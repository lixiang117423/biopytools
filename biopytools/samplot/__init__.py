"""
Samplot Module
Samplot结构变异可视化模块

功能: 使用samplot对基因组结构变异(SV)进行可视化
      支持单区域绘图(plot)和VCF批量绘图(vcf)两种模式

使用示例|Usage Examples:
    from biopytools.samplot import SamplotPlotter, SamplotVcfPlotter

    # 单区域绘图|Single region plot
    plotter = SamplotPlotter(
        bams=["sample.bam"],
        chrom="chr1",
        start=1000,
        end=5000,
        sv_type="DEL",
        output_dir="output"
    )
    plotter.run()

    # VCF批量绘图|Batch plot from VCF
    vcf_plotter = SamplotVcfPlotter(
        bams=["sample.bam"],
        vcf="variants.vcf",
        output_dir="output"
    )
    vcf_plotter.run()

作者|Author: Xiang Li
版本|Version: 1.0.0
日期|Date: 2026-05-11
"""

__version__ = "1.0.0"
__author__ = "Xiang Li"

from .main import SamplotPlotter, SamplotVcfPlotter
from .config import SamplotPlotConfig, SamplotVcfConfig

__all__ = [
    'SamplotPlotter',
    'SamplotVcfPlotter',
    'SamplotPlotConfig',
    'SamplotVcfConfig',
]
