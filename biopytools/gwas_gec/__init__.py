"""
GWAS GEC分析工具包|GWAS GEC Analysis Toolkit
功能: 基于GEC算法计算GWAS显著性阈值，使用VCF格式参考文件|
Features: Calculate GWAS significance thresholds using GEC with VCF format reference file
作者|Author: Xiang LI
版本|Version: 2.2.0
日期|Date: 2026-01-14

使用示例|Usage Examples:
    from biopytools.gwas_gec import GECAnalyzer, GECConfig

    # 使用VCF文件进行GEC分析|Use VCF file for GEC analysis
    analyzer = GECAnalyzer(
        pfile="gwas_results.txt",
        reference="input.vcf.gz",  # VCF格式的参考文件|VCF format reference file
        output_dir="./gec_output",
        threads=12,
        memory="100G"
    )

    # 运行分析|Run analysis
    results = analyzer.run_analysis()
"""

__version__ = "2.3.0"
__author__ = "Xiang LI"

from .main import GECAnalyzer
from .config import GECConfig
from .calculator import GECThresholdCalculator

__all__ = ['GECAnalyzer', 'GECConfig', 'GECThresholdCalculator']
