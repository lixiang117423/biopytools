#!/usr/bin/env python3
"""
VCF连锁不平衡热图生成运行脚本 | VCF LD Heatmap Generator Runner Script
这是一个简化的入口脚本，用于运行VCF连锁不平衡热图分析 | Simple entry script for running VCF LD heatmap analysis

用法 | Usage:
    python run_vcf_ld_heatmap.py -i variants.vcf -o ld_heatmap.png
    
示例 | Examples:
    # 基本LD热图分析 | Basic LD heatmap analysis
    python run_vcf_ld_heatmap.py -i variants.vcf -o ld_heatmap.png
    
    # 指定基因组区域 | Specify genomic region
    python run_vcf_ld_heatmap.py -i variants.vcf -o region_ld.png \\
        --region chr1:1000000-2000000
    
    # 自定义过滤参数 | Custom filtering parameters
    python run_vcf_ld_heatmap.py -i variants.vcf.gz -o ld_map.pdf \\
        --maf 0.05 --max-snps 500 --verbose
    
    # 高质量可视化 | High-quality visualization
    python run_vcf_ld_heatmap.py -i data.vcf -o heatmap.png \\
        --figsize 12 10 --dpi 600 --colormap plasma
    
    # 保存数值矩阵 | Save numeric matrix
    python run_vcf_ld_heatmap.py -i variants.vcf -o heatmap.png \\
        --save-matrix ld_matrix.csv --triangle-only
    
    # 指定样本子集 | Specify sample subset
    python run_vcf_ld_heatmap.py -i variants.vcf -o subset_ld.png \\
        --samples sample1 sample2 sample3 --verbose
    
    # LD阈值过滤 | LD threshold filtering
    python run_vcf_ld_heatmap.py -i variants.vcf -o filtered_ld.png \\
        --ld-threshold 0.2 --title "High LD Regions"
    
    # 完整参数示例 | Complete parameter example
    python run_vcf_ld_heatmap.py -i variants.vcf.gz -o comprehensive_ld.pdf \\
        --region chr2:5000000-6000000 \\
        --maf 0.01 --max-snps 1000 \\
        --figsize 15 12 --dpi 300 \\
        --colormap RdYlBu_r \\
        --title "Chromosome 2 LD Pattern" \\
        --save-matrix chr2_ld_matrix.csv \\
        --ld-threshold 0.1 \\
        --triangle-only --verbose
"""

from biopytools.vcf_ld_heatmap.main import main

if __name__ == "__main__":
    main()
