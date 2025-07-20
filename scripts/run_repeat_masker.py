#!/usr/bin/env python3
"""
基因组重复序列分析运行脚本 | Genome Repeat Sequence Analysis Runner Script
这是一个简化的入口脚本，用于运行基因组重复序列分析 | Simple entry script for genome repeat sequence analysis

用法 | Usage:
    python run_repeat_masker.py -g genome.fasta -s human -o results
    
示例 | Examples:
    # 基本分析 | Basic analysis
    python run_repeat_masker.py -g genome.fa -s human -o repeat_results
    
    # 使用自定义重复库 | Use custom repeat library
    python run_repeat_masker.py -g genome.fa -s mouse -l custom_repeats.lib -o results
    
    # 跳过RepeatModeler，只运行RepeatMasker | Skip RepeatModeler, only run RepeatMasker
    python run_repeat_masker.py -g genome.fa -s human --no-modeler -o results
    
    # 完整分析with更多线程 | Full analysis with more threads
    python run_repeat_masker.py -g genome.fa -s plant -t 16 -m 16 -e 16 -o plant_repeats
    
    # 使用EDTA进行植物基因组分析 | Use EDTA for plant genome analysis
    python run_repeat_masker.py -g maize_genome.fa -s maize --edta-species maize -e 20 -o maize_analysis
    
    # 只运行EDTA分析 | Only run EDTA analysis
    python run_repeat_masker.py -g genome.fa --no-modeler --no-trf --edta-species others -o edta_only
"""

from biopytools.repeat_masker.main import main

if __name__ == "__main__":
    main()
