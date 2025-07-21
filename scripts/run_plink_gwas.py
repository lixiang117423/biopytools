#!/usr/bin/env python3
"""
PLINK GWAS分析运行脚本 | PLINK GWAS Analysis Runner Script
这是一个简化的入口脚本，用于运行PLINK GWAS分析 | Simple entry script for running PLINK GWAS analysis

用法 | Usage:
    python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative -o results
    
示例 | Examples:
    # 质量性状分析 (0/1 -> 1/2转换) | Qualitative trait analysis (0/1 -> 1/2 conversion)
    python run_plink.py -v my_data.vcf.gz -p my_pheno.txt -t qualitative -o plink_results
    
    # 数量性状分析 (保持原值) | Quantitative trait analysis (keep original values)
    python run_plink.py -v my_data.vcf.gz -p my_pheno.txt -t quantitative -o plink_results
    
    # 使用所有校正方法 (默认) | Use all correction methods (default)
    python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative --correction-method all
    
    # 仅使用Bonferroni校正 | Only use Bonferroni correction
    python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative --correction-method bonferroni
    
    # 仅使用FDR校正 | Only use FDR correction
    python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative --correction-method fdr --fdr-alpha 0.1
    
    # 自定义校正参数 | Custom correction parameters
    python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative --correction-method all \\
        --bonferroni-alpha 0.01 --suggestive-threshold 5e-6 --fdr-alpha 0.05
    
    # 自定义质控参数 | Custom QC parameters
    python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative -o results \\
        --mind 0.1 --geno 0.1 --maf 0.05 --hwe 1e-5
    
    # 调整主成分分析参数 | Adjust PCA parameters
    python run_plink.py -v data.vcf.gz -p pheno.txt -t quantitative -o results \\
        --pca-components 15 --pca-use 8
    
    # 使用多线程 | Use multiple threads
    python run_plink.py -v data.vcf.gz -p pheno.txt -t qualitative -o results --threads 8
"""

from biopytools.plink_gwas.main import main

if __name__ == "__main__":
    main()
