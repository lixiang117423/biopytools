#!/usr/bin/env python3
"""
K-mer数据库分析运行脚本 | K-mer Database Analysis Runner Script
这是一个简化的入口脚本，用于运行k-mer数据库分析 | Simple entry script for running k-mer database analysis

用法 | Usage:
    python run_kmer.py -g genes.fasta -f /path/to/fastq -o ./results
    
示例 | Examples:
    # 首次运行 (构建数据库) | First run (build database)
    python run_kmer.py -g genes.fasta -f /path/to/fastq -o ./db_results -k 51 -t 32
    
    # 后续查询新基因 (跳过构建) | Subsequent queries of new genes (skip build)
    python run_kmer.py -g new_genes.fasta -f /path/to/fastq -o ./db_results --skip-build
    
    # 运行单倍型分析 | Run haplotype analysis
    python run_kmer.py -g genes.fasta -f /path/to/fastq -o ./results --run-haplotype -n 5
    
    # 大规模数据集优化 | Large dataset optimization
    python run_kmer.py -g genes.fasta -f /path/to/fastq -o ./results -t 64 -m 3
    
    # 调试模式 | Debug mode
    python run_kmer.py -g genes.fasta -f /path/to/fastq -o ./results --debug
"""

from biopytools.kmer_old.main import main

if __name__ == "__main__":
    main()
