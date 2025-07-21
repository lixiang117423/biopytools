#!/usr/bin/env python3
"""
序列提取运行脚本 | Sequence Extraction Runner Script
这是一个简化的入口脚本，用于运行序列提取分析 | Simple entry script for running sequence extraction analysis

用法 | Usage:
    python run_sequence_extractor.py -v variants.vcf -g genome.fa -c chr1 -s 1000 -e 1050
    
示例 | Examples:
    # 基本提取 | Basic extraction
    python run_sequence_extractor.py -v variants.vcf -g genome.fa -c chr1 -s 1000 -e 1050
    
    # 指定输出格式和目录 | Specify output format and directory
    python run_sequence_extractor.py -v variants.vcf.gz -g genome.fa -c chr1 -s 1000 -e 1050 \\
        -o results --format fasta
    
    # 使用第二等位基因并排除特定样品 | Use second allele and exclude specific samples
    python run_sequence_extractor.py -v variants.vcf -g genome.fa -c chr1 -s 1000 -e 1050 \\
        --second-allele --exclude-samples "sample1,sample2"
    
    # 质量过滤并指定样品 | Quality filtering and specify samples
    python run_sequence_extractor.py -v variants.vcf -g genome.fa -c chr1 -s 1000 -e 1050 \\
        --min-qual 30 --samples samples.txt
"""

from biopytools.vcf_sequence.main import main

if __name__ == "__main__":
    main()
