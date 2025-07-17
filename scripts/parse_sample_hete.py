#!/usr/bin/env python3
"""
VCF基因型统计运行脚本 | VCF Genotype Statistics Runner Script
这是一个简化的入口脚本，用于运行VCF基因型统计分析 | Simple entry script for running VCF genotype statistics analysis

支持的基因型格式 | Supported genotype formats:
- 未定相: 0/0, 0/1, 1/1, ./. | Unphased: 0/0, 0/1, 1/1, ./.
- 已定相: 0|0, 0|1, 1|1, .|. | Phased: 0|0, 0|1, 1|1, .|.
- 多等位基因: 0/2, 1/2等 | Multi-allelic: 0/2, 1/2, etc.

用法 | Usage:
    python run_vcf_stats.py -v variants.vcf -o vcf_stats_output
    
示例 | Examples:
    # 基本分析 | Basic analysis
    python run_vcf_stats.py -v my_variants.vcf -o my_stats
    
    # 应用过滤条件 | Apply filtering conditions
    python run_vcf_stats.py -v variants.vcf.gz -o filtered_stats -d 10 -q 30.0 -e
    
    # 仅输出汇总统计 | Summary statistics only
    python run_vcf_stats.py -v variants.vcf -o simple_stats -D
    
    # 长参数格式 | Long parameter format
    python run_vcf_stats.py --vcf variants.vcf --output filtered_stats \\
        --min-depth 10 --min-qual 30.0 --exclude-missing
"""

from biopytools.vcf_stats_sample.main import main

if __name__ == "__main__":
    main()
