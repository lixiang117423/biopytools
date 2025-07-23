"""
PLINK GWAS分析工具包 | PLINK GWAS Analysis Toolkit
功能: VCF到GWAS分析的完整流程，支持二分类性状和群体分层控制 | 
Features: Complete pipeline from VCF to GWAS analysis, supporting binary traits and population stratification
作者 | Author: Xiang LI  
版本 | Version: v10 - 模块化重构版 | Modular refactored version
日期 | Date: 2025-07-14

使用示例 | Usage Examples:
    from plink import PlinkGWAS, PlinkGWASConfig
    
    # 创建分析器 | Create analyzer
    gwas = PlinkGWAS(
        vcf_file="data.vcf.gz",
        phenotype_file="pheno.txt",
        output_dir="gwas_results"
    )
    
    # 运行分析 | Run analysis
    gwas.run_analysis()
"""

__version__ = "10.0.0"
__author__ = "Xiang LI"

from .main import PlinkGWAS
from .config import PlinkGWASConfig

__all__ = ['PlinkGWAS', 'PlinkGWASConfig']
