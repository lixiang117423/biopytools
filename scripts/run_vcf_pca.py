#!/usr/bin/env python3
"""
VCF PCA分析运行脚本 | VCF PCA Analysis Runner Script
这是一个简化的入口脚本，用于运行VCF PCA分析 | Simple entry script for running VCF PCA analysis

用法 | Usage:
    python run_vcf_pca.py --vcf variants.vcf --output pca_results
    
示例 | Examples:
    # 基本PCA分析 | Basic PCA analysis
    python run_vcf_pca.py -v variants.vcf -o pca_results
    
    # 指定主成分数量和可视化 | Specify PC number and visualization
    python run_vcf_pca.py -v data.vcf.gz -o results \\
        -c 15 -p
    
    # 包含样本信息和分组 | Include sample info and grouping
    python run_vcf_pca.py -v variants.vcf -o pca_out \\
        -s samples.txt -g population -p
    
    # 跳过质控过滤 | Skip quality control filtering
    python run_vcf_pca.py -v filtered_variants.vcf -o results \\
        --skip-qc -p
    
    # 自定义质控参数 | Custom QC parameters
    python run_vcf_pca.py -v data.vcf -o results \\
        -m 0.01 --missing 0.05 -c 20
"""

from biopytools.vcf_pca.main import main

if __name__ == "__main__":
    main()
