"""
PLINK GWAS分析工具包 | PLINK GWAS Analysis Toolkit
功能: VCF到GWAS分析的完整流程，支持质量性状和数量性状分析，多种显著性校正方法 | 
Features: Complete pipeline from VCF to GWAS analysis, supporting both qualitative and quantitative traits, multiple significance correction methods
作者 | Author: Xiang LI  
版本 | Version: v12 - 支持多种显著性校正方法 | Supporting multiple significance correction methods
日期 | Date: 2025-07-16

使用示例 | Usage Examples:
    from plink import PlinkGWAS, PlinkGWASConfig
    
    # 质量性状分析，使用所有校正方法 | Qualitative trait analysis with all correction methods
    gwas = PlinkGWAS(
        vcf_file="data.vcf.gz",
        phenotype_file="pheno.txt",
        trait_type="qualitative",
        correction_method="all",
        output_dir="gwas_results"
    )
    
    # 数量性状分析，仅使用FDR校正 | Quantitative trait analysis with FDR correction only
    gwas = PlinkGWAS(
        vcf_file="data.vcf.gz",
        phenotype_file="pheno.txt",
        trait_type="quantitative",
        correction_method="fdr",
        fdr_alpha=0.05,
        output_dir="gwas_results"
    )
    
    # 运行分析 | Run analysis
    gwas.run_analysis()
"""

__version__ = "12.0.0"
__author__ = "Xiang LI"

from .main import PlinkGWAS
from .config import PlinkGWASConfig

__all__ = ['PlinkGWAS', 'PlinkGWASConfig']
