#!/usr/bin/env python3
"""
VCF系统发育分析运行脚本 v2.0 | VCF Phylogenetic Analysis Runner Script v2.0
这是一个简化的入口脚本，用于运行VCF系统发育分析 | Simple entry script for running VCF phylogenetic analysis

用法 | Usage:
    python run_vcf_phylo.py --input variants.vcf --output phylo_results
    
示例 | Examples:
    # 基本系统发育分析 | Basic phylogenetic analysis
    python run_vcf_phylo.py -i variants.vcf -o phylo_results
    
    # 指定完整路径 | Specify full paths
    python run_vcf_phylo.py -i 01.data/wild.snp.new.record.vcf.recode.vcf \\
        -o 02.tree/wild_snp_dis_mat -t 02.tree/wild_snp_dis.nwk
    
    # 从已有距离矩阵构建树 | Build tree from existing matrix
    python run_vcf_phylo.py -d existing_matrix.txt -t tree.nwk --skip-vcf2dis
    
    # 自定义VCF2Dis路径 | Custom VCF2Dis path
    python run_vcf_phylo.py -i data.vcf -o results --vcf2dis-path /path/to/VCF2Dis

更新日志 | Changelog:
    v2.0: 简化实现，直接使用scikit-bio构建NJ树
    v1.0: 初始版本，支持多种树构建方法
"""

from biopytools.vcf_phylo.main import main

if __name__ == "__main__":
    main()
