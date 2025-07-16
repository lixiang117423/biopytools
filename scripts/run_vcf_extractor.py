#!/usr/bin/env python3
"""
VCF基因型提取运行脚本 | VCF Genotype Extraction Runner Script
这是一个简化的入口脚本，用于运行VCF基因型提取 | Simple entry script for running VCF genotype extraction

用法 | Usage:
    python run_vcf_extractor.py variants.vcf.gz --output genotype_results
    
示例 | Examples:
    # 基本提取（所有样本，TXT格式） | Basic extraction (all samples, TXT format)
    python run_vcf_extractor.py variants.vcf.gz -o genotype_results
    
    # 指定特定样本 | Extract specific samples
    python run_vcf_extractor.py variants.vcf.gz -s "OV8-1,OV8-105,OV8-106" -o selected_samples
    
    # 只保留双等位位点，按染色体拆分输出 | Only biallelic sites, split by chromosome
    python run_vcf_extractor.py variants.vcf.gz --biallelic-only --each -o biallelic_by_chrom
    
    # 输出CSV格式 | Output CSV format
    python run_vcf_extractor.py variants.vcf.gz -t csv -o genotype_results
    
    # 输出Excel格式（需要pandas） | Output Excel format (requires pandas)
    python run_vcf_extractor.py variants.vcf.gz -t excel -o genotype_results
    
    # 指定输出目录 | Specify output directory
    python run_vcf_extractor.py variants.vcf.gz -o results --output-dir ./output_folder
"""

from biopytools.vcf_genotype_extractor.main import main

if __name__ == "__main__":
    main()
