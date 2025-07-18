#!/usr/bin/env python3
"""
K-mer PAV分析运行脚本 (双阶段设计，含from列) | K-mer PAV Analysis Runner Script (Two-Stage Design with from column)
这是一个简化的入口脚本，用于运行K-mer存在/缺失变异分析 | 
Simple entry script for running K-mer Presence/Absence Variation analysis

双阶段设计 | Two-Stage Design:
  阶段1 | Phase 1: 从数据库文件构建统一k-mer数据库，同时追踪k-mer来源
  阶段2 | Phase 2: 查询文件与k-mer数据库比较分析

样本处理逻辑 | Sample Processing Logic:
  - FASTQ文件: 整个文件作为一个样本
  - FASTA文件: 每条序列作为一个样本

新增功能 | New Features:
  - from列: 显示k-mer来源于database中的哪个文件
  - feature列: 显示k-mer在query样本中的分布模式

依赖工具 | Required Tools:
    - KMC: 高效的k-mer计数工具 | Efficient k-mer counting tool
    - kmc_tools: KMC工具套件 | KMC tools suite
    - BioPython: FASTA序列处理 | FASTA sequence processing

安装方法 | Installation:
    Ubuntu/Debian: sudo apt-get install kmc
    Conda: conda install -c bioconda kmc
    Python: pip install biopython pandas numpy

用法 | Usage:
    python run_kmer_pav.py --database-input db_files/ --query-input samples/ -o analysis
    
示例 | Examples:
    # 基本分析 | Basic analysis
    python run_kmer_pav.py --database-input database_files/ --query-input query_samples/ -o my_analysis
    
    # 单个数据库文件，多个查询文件 | Single database file, multiple query files
    python run_kmer_pav.py --database-input reference.fasta --query-input samples/ -o analysis
    
    # 多个数据库文件，单个查询文件 | Multiple database files, single query file
    python run_kmer_pav.py --database-input db_files/ --query-input query.fasta -s 25
    
    # 指定文件模式 | Specify file patterns
    python run_kmer_pav.py \\
        --database-input db_files/ --database-pattern "*.reference.fa" \\
        --query-input samples/ --query-pattern "*_clean.fq.gz" \\
        -o analysis
    
    # 完整参数示例 | Full parameter example
    python run_kmer_pav.py \\
        --database-input database_files/ \\
        --query-input query_samples/ \\
        -o kmer_analysis \\
        -s 31 -r \\
        --min-count 2 --max-count 10000 \\
        --threads 16 --kmc-memory 32 \\
        --output-dir results/ --keep-intermediate

输出文件说明 | Output File Description:
    主要输出文件包含两个重要列:
    - from列: 显示k-mer在database文件中的来源
      * common: 所有database文件都包含
      * 文件名: 仅该database文件包含  
      * variable(x/y): 在x个database文件中包含(共y个文件)
    - feature列: 显示k-mer在query样本中的分布
      * common: 所有query样本都包含
      * 样本名: 仅该query样本包含
      * variable(x/y): 在x个query样本中包含(共y个样本)
"""

from biopytools.kmer_pav.main import main

if __name__ == "__main__":
    main()
