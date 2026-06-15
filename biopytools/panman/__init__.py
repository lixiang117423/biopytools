"""
PanMAN泛基因组分析工具包|PanMAN Pangenome Analysis Toolkit
功能: 构建和分析PanMAN（泛基因组突变注释网络）|Features: Build and analyze PanMANs (Pangenome Mutation-Annotated Networks)
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-01-15

使用示例|Usage Examples:
    from biopytools.panman import PanMANBuildRunner, PanMANExtractRunner, PanMANConfig

    # 构建PanMAN|Build PanMAN
    builder = PanMANBuildRunner(
        pangraph_file="input.json",
        newick_file="tree.nwk",
        output_prefix="my_panman"
    )
    builder.run_analysis()

    # 提取数据|Extract data
    extractor = PanMANExtractRunner(
        panman_file="data.panman",
        output_prefix="output",
        extract_summary=True,
        extract_fasta=True,
        extract_vcf=True
    )
    extractor.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import PanMANBuildRunner, PanMANExtractRunner
from .config import PanMANConfig

__all__ = [
    'PanMANBuildRunner',
    'PanMANExtractRunner',
    'PanMANConfig'
]
