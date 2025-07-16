#!/usr/bin/env python3
"""
GFF3基因转录本提取运行脚本 | GFF3 Gene Transcript Extraction Runner Script
这是一个简化的入口脚本，用于运行GFF3基因转录本提取分析 | Simple entry script for running GFF3 gene transcript extraction analysis

用法 | Usage:
    python run_gff_extractor.py -g annotation.gff3 -o gene_transcript_info.tsv
    
示例 | Examples:
    # 基本提取 | Basic extraction
    python run_gff_extractor.py -g input.gff3 -o output.tsv
    
    # 指定基因和转录本类型 | Specify gene and transcript types
    python run_gff_extractor.py -g input.gff3 -o output.tsv --gene-type gene --transcript-types mRNA transcript
    
    # 处理特定类型的转录本 | Process specific transcript types
    python run_gff_extractor.py -g input.gff3 -o output.tsv --transcript-types mRNA
"""

from biopytools.gff_utils.main import main

if __name__ == "__main__":
    main()
