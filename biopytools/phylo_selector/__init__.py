"""
系统发育树样品选择工具包|Phylogenetic Tree Sample Selector Toolkit
功能: 从系统发育树中选择代表性样品，支持分层抽样和最大间距优化|
Features: Select representative samples from phylogenetic trees, supporting stratified sampling and maximum spacing optimization
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-01-30

使用示例|Usage Examples:
    from biopytools.phylo_selector import PhyloSelectorRunner, PhyloSelectorConfig

    # 简单选择（不分组）
    runner = PhyloSelectorRunner(
        newick_file="tree.nwk",
        output_prefix="selected_samples"
    )
    runner.run()

    # 分层选择（带分组信息）
    runner = PhyloSelectorRunner(
        newick_file="tree.nwk",
        group_file="groups.txt",
        output_prefix="selected_samples",
        n_samples=100
    )
    runner.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import PhyloSelectorRunner
from .config import PhyloSelectorConfig
from .calculator import PhyloSelectorCalculator
from .parser import NewickParser, GroupTableParser, PCAFileParser
from .utils import PhyloSelectorLogger, format_number

__all__ = [
    'PhyloSelectorRunner',
    'PhyloSelectorConfig',
    'PhyloSelectorCalculator',
    'NewickParser',
    'GroupTableParser',
    'PCAFileParser',
    'PhyloSelectorLogger',
    'format_number'
]
