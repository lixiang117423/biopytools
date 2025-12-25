"""
🌳 IQ-TREE 系统发育树构建工具包 | IQ-TREE Phylogenetic Tree Construction Toolkit
功能: 基于最大似然法的系统发育树构建，支持模型选择、Bootstrap分析和多种高级功能 |
Features: Maximum likelihood phylogenetic tree construction with model selection, bootstrap analysis and advanced features
作者 | Author: Claude
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-10-06

使用示例 | Usage Examples:
    from biopytools.iqtree import IQTreeAnalyzer, TreeConfig
    
    # 创建分析器 | Create analyzer
    analyzer = IQTreeAnalyzer(
        input_file="alignment.fasta",
        output_dir="tree_results",
        prefix="my_tree",
        bootstrap=1000,
        threads=88
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import IQTreeAnalyzer
from .config import TreeConfig

__all__ = ['IQTreeAnalyzer', 'TreeConfig']
