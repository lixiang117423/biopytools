#!/usr/bin/env python3
"""
最长转录本提取运行脚本 | Longest mRNA Extraction Runner Script
这是一个简化的入口脚本，用于运行最长转录本提取 | Simple entry script for running longest mRNA extraction

用法 | Usage:
    python run_longest_mrna.py -g genome.fa -f annotation.gff3 -o longest_transcripts.fa
    
示例 | Examples:
    # 基本提取 | Basic extraction
    python run_longest_mrna.py -g genome.fa -f genes.gff3 -o longest_transcripts.fa
    
    # 指定基因信息输出文件 | Specify gene info output file
    python run_longest_mrna.py -g genome.fa -f genes.gff3 -o longest_transcripts.fa --gene-info gene_info.txt
"""

from biopytools.longest_mrna.main import main

if __name__ == "__main__":
    main()
