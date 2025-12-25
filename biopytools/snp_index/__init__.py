"""
SNP index计算和分析工具包 | SNP Index Calculation and Analysis Toolkit
功能: 从VCF文件计算SNP index和ΔSNP index，提供结果分析和可视化功能
Features: Calculate SNP index and ΔSNP index from VCF files, provide result analysis and visualization

作者 | Author: Xiang LI
版本 | Version: 1.0.0
日期 | Date: 2025-12-20

使用示例 | Usage Examples:
    from biopytools.snp_index import SNPIndexCalculator, SNPIndexAnalyzer

    # 计算SNP index | Calculate SNP index
    calculator = SNPIndexCalculator(
        input_vcf="variants.vcf.gz",
        output_file="results.tsv",
        min_depth=10,
        min_quality=20
    )
    calculator.calculate()

    # 分析结果 | Analyze results
    analyzer = SNPIndexAnalyzer(result_file="results.tsv")
    analyzer.analyze()
    analyzer.visualize(output_prefix="snp_analysis")
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import SNPIndexProcessor
from .config import SNPIndexConfig
from .calculator import SNPIndexCalculator
from .analyzer import SNPIndexAnalyzer
from .visualizer import SNPIndexVisualizer

__all__ = [
    'SNPIndexProcessor',
    'SNPIndexConfig',
    'SNPIndexCalculator',
    'SNPIndexAnalyzer',
    'SNPIndexVisualizer'
]