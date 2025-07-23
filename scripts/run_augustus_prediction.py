#!/usr/bin/env python3
"""
Augustus基因预测流水线运行脚本 | Augustus Gene Prediction Pipeline Runner Script
这是一个简化的入口脚本，用于运行Augustus基因预测流水线 | Simple entry script for running Augustus gene prediction pipeline

用法 | Usage:
    python run_augustus_pipeline.py -s Rice_NLR_Model -g genome.fa -a annotations.gff3
    
示例 | Examples:
    # 基本用法 | Basic usage
    python run_augustus_pipeline.py -s Rice_NLR_Model -g genome.fa -a annotations.gff3
    
    # 自定义输出目录和训练比例 | Custom output directory and training ratio
    python run_augustus_pipeline.py -s MySpecies -g genome.fa -a genes.gff3 -o my_results -t 0.85
    
    # 指定Augustus路径和侧翼长度 | Specify Augustus path and flank length
    python run_augustus_pipeline.py -s Plant_Model -g plant.fa -a plant.gff3 -p /opt/augustus/bin -f 1500
    
    # 完整参数示例 | Complete parameter example
    python run_augustus_pipeline.py \\
        --species-name Rice_NLR_Model \\
        --genome-file genome.fa \\
        --gff-file annotations.gff3 \\
        --output-dir ./augustus_results \\
        --train-ratio 0.8 \\
        --flank-length 1000 \\
        --augustus-path /path/to/augustus/bin
"""

from biopytools.augustus_prediction.main import main

if __name__ == "__main__":
    main()
