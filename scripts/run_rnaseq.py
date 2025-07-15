#!/usr/bin/env python3
"""
RNA-seq分析运行脚本 | RNA-seq Analysis Runner Script
这是一个简化的入口脚本，用于运行RNA-seq分析 | Simple entry script for running RNA-seq analysis

用法 | Usage:
    python run_rnaseq.py -g genome.fa -f annotation.gtf -i fastq_dir -o results
    
示例 | Examples:
    # 基本分析 | Basic analysis
    python run_rnaseq.py -g genome.fa -f genes.gtf -i ./fastq -o rnaseq_results
    
    # 指定文件模式和线程数 | Specify file pattern and thread count
    python run_rnaseq.py -g genome.fa -f genes.gtf -i ./fastq -o results \\
        -p "*.R1.fastq.gz" -t 16
    
    # 删除BAM文件以节省空间 | Remove BAM files to save space
    python run_rnaseq.py -g genome.fa -f genes.gtf -i ./fastq -o results -r yes
    
    # 使用样本信息文件 | Use sample information file
    python run_rnaseq.py -g genome.fa -f genes.gtf -i samples.txt -o results
"""

from biopytools.rnaseq.main import main

if __name__ == "__main__":
    main()
