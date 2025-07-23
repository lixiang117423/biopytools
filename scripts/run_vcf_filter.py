#!/usr/bin/env python3
"""
VCF文件筛选运行脚本 | VCF File Filtering Runner Script
这是一个简化的入口脚本，用于运行VCF文件筛选 | Simple entry script for running VCF file filtering

用法 | Usage:
    python run_vcf_filter.py -i variants.vcf -c chr1 -s 1000 -e 2000
    
示例 | Examples:
    # 基本筛选（高性能模式，默认跳过验证）| Basic filtering (high performance mode, skip validation by default)
    python run_vcf_filter.py -i variants.vcf -c chr1 -s 1000 -e 2000
    
    # 多染色体筛选 | Multiple chromosome filtering
    python run_vcf_filter.py -i variants.vcf -c "chr1,chr2,chr3" -o multi_chr.vcf
    
    # 使用PLINK格式转换 | Using PLINK format conversion
    python run_vcf_filter.py -i variants.vcf -c chr1 --convert-format --maf 0.05 --max-missing 0.1
    
    # 质量筛选 | Quality filtering
    python run_vcf_filter.py -i variants.vcf -c chr1 --quality-threshold 30 --biallelic-only
    
    # 样本筛选 | Sample filtering
    python run_vcf_filter.py -i variants.vcf -c chr1 --keep-samples "sample1,sample2,sample3"
    
    # 变异类型筛选 | Variant type filtering
    python run_vcf_filter.py -i variants.vcf -c chr1 --biallelic-only --remove-indels
    
    # 综合筛选（高性能）| Comprehensive filtering (high performance)
    python run_vcf_filter.py -i variants.vcf.gz -c chr1 -s 1000000 -e 5000000 \\
        --quality-threshold 30 --maf 0.01 \\
        --biallelic-only --remove-indels \\
        --keep-samples "sample1,sample2,sample3" \\
        --verbose
    
    # PLINK高级筛选 | PLINK advanced filtering
    python run_vcf_filter.py -i variants.vcf.gz -c chr1 \\
        --convert-format --plink-path /usr/local/bin/plink \\
        --maf 0.05 --max-missing 0.05 --allow-extra-chr \\
        --verbose
        
    # 如果需要验证（较慢但更安全）| If validation needed (slower but safer)
    python run_vcf_filter.py -i variants.vcf -c chr1 --force-validation --verbose
"""

from biopytools.vcf_filter.main import main

if __name__ == "__main__":
    main()
