#!/usr/bin/env python3
"""
FASTP质控运行脚本 | FASTP Quality Control Runner Script
这是一个简化的入口脚本，用于运行FASTQ质控分析 | Simple entry script for running FASTQ quality control analysis

用法 | Usage:
    python run_fastp.py -i raw_data -o clean_data
    
示例 | Examples:
    # 基本质控 | Basic quality control
    python run_fastp.py -i ./raw_data -o ./clean_data
    
    # 自定义参数 | Custom parameters
    python run_fastp.py -i ./raw_data -o ./clean_data \\
        --fastp-path /path/to/fastp --threads 16 --quality-threshold 25
    
    # 指定文件后缀 | Specify file suffixes
    python run_fastp.py -i ./raw_data -o ./clean_data \\
        --read1-suffix "_R1.fastq.gz" --read2-suffix "_R2.fastq.gz"
"""

from biopytools.fastp.main import main

if __name__ == "__main__":
    main()
