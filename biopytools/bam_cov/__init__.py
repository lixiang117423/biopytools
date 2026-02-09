"""
BAM覆盖度统计工具包|BAM Coverage Statistics Toolkit
功能: 计算BAM文件中指定区域的碱基覆盖度，支持单个或多个BAM文件处理
Features: Calculate per-base coverage in specified regions from BAM files, support single or multiple BAM files
作者|Author: Xiang LI
版本|Version: v1.0.0
日期|Date: 2025-12-31

使用示例|Usage Examples:
    from biopytools.bam_coverage_stats import BAMCoverageAnalyzer, BAMCoverageConfig

    # 创建分析器|Create analyzer
    analyzer = BAMCoverageAnalyzer(
        bam_dir="./bam_files",
        chromosome="chr1",
        start=1000000,
        end=2000000,
        output_prefix="coverage_stats"
    )

    # 运行分析|Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import BAMCoverageAnalyzer
from .config import BAMCoverageConfig

__all__ = ['BAMCoverageAnalyzer', 'BAMCoverageConfig']
