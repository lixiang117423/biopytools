#!/usr/bin/env python3
"""
PacBio HiFi结构变异检测运行脚本 | PacBio HiFi Structural Variant Detection Runner Script
这是一个简化的入口脚本，用于运行PacBio HiFi SV检测分析 | Simple entry script for running PacBio HiFi SV detection analysis

用法 | Usage:
    python run_pacbio_sv.py -b sample.bam -r reference.fa -s sample_name
    
示例 | Examples:
    # 基本分析 | Basic analysis
    python run_pacbio_sv.py -b OV53-1.bam -r ref_genome.fa -s OV53-1
    
    # 指定输出目录和线程数 | Specify output directory and thread count
    python run_pacbio_sv.py -b sample.bam -r genome.fa -s sample01 \\
        -o sv_results -t 32
    
    # 自定义SV检测参数 | Custom SV detection parameters
    python run_pacbio_sv.py -b data.bam -r ref.fa -s test \\
        --min-sv-length 100 --min-support 5 -q 30
    
    # 大片段SV分析 | Large SV analysis
    python run_pacbio_sv.py -b aligned.bam -r reference.fa -s sample \\
        --large-sv-threshold 5000 --very-large-sv-threshold 50000
    
    # 完整参数示例 | Complete parameter example
    python run_pacbio_sv.py -b sample.bam -r genome.fa -s test_sample \\
        -o custom_output -t 24 --min-sv-length 75 --min-support 4 \\
        -q 25 --large-sv-threshold 2000 --survivor-distance 500
"""

from biopytools.pacbio_sv.main import main

if __name__ == "__main__":
    main()
