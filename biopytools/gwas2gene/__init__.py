"""
GWAS2Gene - GWAS候选基因筛选工具|GWAS Candidate Gene Finder
功能: 根据GWAS结果筛选显著SNP，在指定窗口内从GFF文件中提取候选基因|
Features: Extract candidate genes near significant GWAS SNPs from GFF annotations
作者|Author: Xiang LI
版本|Version: v1.0 - 初始版本|Initial version
日期|Date: 2026-03-23

使用示例|Usage Examples:
    from biopytools.gwas2gene import GWAS2GeneFinder, GWAS2GeneConfig

    # 基本用法|Basic usage
    finder = GWAS2GeneFinder(
        gwas_file='gwas_result.txt',
        pval_col='Pvalue',
        threshold=1e-5,
        window=200000,
        gff_file='annotation.gff3',
        output_file='candidate_genes.tsv'
    )

    finder.run()

    # 带功能注释|With function annotations
    finder = GWAS2GeneFinder(
        gwas_file='gwas_result.txt',
        pval_col='Pvalue',
        threshold=1e-5,
        window=200000,
        gff_file='annotation.gff3',
        output_file='candidate_genes.tsv',
        func_file='gene_function.txt'
    )

    finder.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import GWAS2GeneFinder
from .config import GWAS2GeneConfig

__all__ = ['GWAS2GeneFinder', 'GWAS2GeneConfig']
