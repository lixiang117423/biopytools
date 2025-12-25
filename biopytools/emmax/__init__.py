"""
🧬 AutoGWAS: 自动化GWAS分析工具包 | AutoGWAS: Automated GWAS Analysis Toolkit
功能: 基于EMMAX的完整GWAS分析流程，支持曼哈顿图、QQ图生成 |
Features: Complete GWAS analysis pipeline based on EMMAX, supporting Manhattan plots and QQ plots
作者 | Author: MiniMax Agent
版本 | Version: v1.0.0
日期 | Date: 2025-11-04

使用示例 | Usage Examples:
    from gwas_analysis import GWASAnalyzer, GWASConfig

    # 创建分析器 | Create analyzer
    analyzer = GWASAnalyzer(
        vcf_file="variants.vcf",
        phenotype_file="phenotypes.txt",
        output_prefix="gwas_results"
    )

    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "MiniMax Agent"

# 条件导入以支持单文件运行
try:
    from .config import GWASConfig
    from .main import GWASAnalyzer
except (ImportError, SystemError):
    # 单文件运行或模块不存在时，跳过导入
    # 这些将在需要时被动态导入
    GWASAnalyzer = None
    GWASConfig = None

__all__ = ["GWASAnalyzer", "GWASConfig"]
