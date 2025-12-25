"""
甲烷循环相关基因丰度分析工具 | Methane Cycle Gene Abundance Analysis Tool
功能: 基于MCycDB数据库分析甲烷循环相关基因的丰度 |
Features: Analysis of methane cycle-related genes based on MCycDB database
作者 | Author: Xiang LI
版本 | Version: 1.0.0 - 模块化版本 | Modular version
日期 | Date: 2025-12-19

使用示例 | Usage Examples:
    from biopytools.mcyc import MCycAnalyzer, MCycConfig

    # 创建分析器 | Create analyzer
    config = MCycConfig(
        input_list="samples.txt",
        output_dir="./mcyc_results"
    )
    analyzer = MCycAnalyzer(config)
    analyzer.run_full_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import MCycAnalyzer
from .config import MCycConfig
from .calculator import MatrixCalculator

__all__ = ['MCycAnalyzer', 'MCycConfig', 'MatrixCalculator']