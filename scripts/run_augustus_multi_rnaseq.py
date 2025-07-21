#!/usr/bin/env python3
"""
多转录组Augustus基因预测运行脚本 | Multiple RNA-seq Augustus Gene Prediction Runner Script
完整修复版本 | Complete Fixed Version

用法 | Usage:
    # 方式1: 自动发现样本 (推荐) | Method 1: Auto-discover samples (Recommended)
    python run_augustus_multi_rnaseq.py -g genome.fasta -i /path/to/fastq -s model_name
    
    # 方式2: 使用配置文件 | Method 2: Using config file
    python run_augustus_multi_rnaseq.py -g genome.fasta -c samples.txt -s model_name
    
示例 | Examples:
    # 基本预测 (自动发现样本) | Basic prediction (auto-discover samples)
    python run_augustus_multi_rnaseq.py -g genome.fasta -i ./fastq_data -s arabidopsis
    
    # 指定文件模式 | Specify file pattern
    python run_augustus_multi_rnaseq.py -g genome.fasta -i ./fastq_data -s arabidopsis \\
        -p "*.R1.fq.gz"
    
    # 完整参数示例 | Full parameter example
    python run_augustus_multi_rnaseq.py -g genome.fasta -i ./fastq_data -s arabidopsis \\
        -o augustus_results -t 16 -p "*.read1.fastq.gz"
    
    # 跳过依赖检查 | Skip dependency check
    python run_augustus_multi_rnaseq.py -g genome.fasta -i ./fastq_data -s arabidopsis \\
        --skip-deps-check

参数说明 | Parameter Description:
    必需参数 | Required:
        -g, --genome        基因组fasta文件 | Genome fasta file
        -s, --species       Augustus物种模型 | Augustus species model
        
    输入方式 (二选一) | Input method (choose one):
        -i, --input-dir     FASTQ文件目录 | FASTQ files directory  
        -c, --config        样本配置文件 | Sample configuration file
    
    可选参数 | Optional:
        -p, --pattern       R1文件匹配模式 | R1 file pattern (default: *.R1.fastq.gz)
        -o, --output        输出目录 | Output directory (default: augustus_multi_rnaseq)
        -t, --threads       线程数 | Thread count (default: 8)
        -x, --hisat2-index  HISAT2索引前缀 | HISAT2 index prefix
        -m, --min-intron-support  最小内含子支持度 | Min intron support (default: 2)
        -f, --no-filter-bam       跳过BAM过滤 | Skip BAM filtering
        -a, --no-alternatives     禁用可变剪切 | Disable alternative splicing
        -z, --splicesites   剪切位点类型 | Splice site types (atac/gtag/gcag, default: atac)
        --skip-deps-check   跳过依赖检查 | Skip dependency check

配置文件格式 | Configuration file format (when using -c):
    sample1_name /path/to/sample1_R1.fastq /path/to/sample1_R2.fastq
    sample2_name /path/to/sample2_R1.fastq /path/to/sample2_R2.fastq
    ...
"""

from biopytools.augustus_multi_rnaseq.main import main

if __name__ == "__main__":
    main()
