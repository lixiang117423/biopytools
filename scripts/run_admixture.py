#!/usr/bin/env python3
"""
ADMIXTURE群体结构分析运行脚本 | ADMIXTURE Population Structure Analysis Runner Script
这是一个简化的入口脚本，用于运行ADMIXTURE分析 | Simple entry script for running ADMIXTURE analysis

用法 | Usage:
    # 使用长参数
    python run_admixture.py --vcf data.vcf.gz --output admixture_results
    
    # 使用短参数
    python run_admixture.py -v data.vcf.gz -o admixture_results -t 8
    
示例 | Examples:
    # 基本分析 | Basic analysis
    python run_admixture.py -v final_filtered.vcf.gz -o results
    
    # 指定K值范围和线程数 | Specify K range and thread count
    python run_admixture.py -v data.vcf.gz -o results \\
        --min-k 2 --max-k 8 -t 8
    
    # 跳过预处理步骤 | Skip preprocessing steps
    python run_admixture.py -v clean_data.vcf.gz -o results \\
        --skip-preprocessing
    
    # 自定义质控参数 | Custom QC parameters
    python run_admixture.py -v data.vcf.gz -o results \\
        --maf 0.05 --missing 0.05 --hwe 1e-5
"""

from biopytools.admixture.main import main

if __name__ == "__main__":
    main()
