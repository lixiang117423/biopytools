#!/usr/bin/env python3
"""
🧬 甲基化分析运行脚本 | Methylation Analysis Runner Script
这是一个简化的入口脚本，用于运行甲基化分析流程

用法 | Usage:
    python run_methylation_pipeline.py -r /path/to/raw -g /path/to/genome.fa -o /path/to/output
    
示例 | Examples:
    # 基础分析 | Basic analysis
    python run_methylation_pipeline.py -r ./raw_data -g ./genome.fa -o ./results
    
    # 增强分析 | Enhanced analysis
    python run_methylation_pipeline.py -r ./raw_data -g ./genome.fa -o ./results \\
        -p ./target_promoter.fa -a ./genome.gff3 -e
    
    # 自定义参数 | Custom parameters
    python run_methylation_pipeline.py -r ./raw_data -g ./genome.fa -o ./results \\
        -j 16 -c 10 -d 0.25 -v 0.001 -e
    
    # 试运行模式（仅检查配置）| Dry run mode (check config only)
    python run_methylation_pipeline.py -r ./raw_data -g ./genome.fa -o ./results \\
        -e --dry-run
    
    # 详细输出模式 | Verbose output mode
    python run_methylation_pipeline.py -r ./raw_data -g ./genome.fa -o ./results \\
        -e --verbose
    
    # 使用自定义样品分组 | Custom sample grouping
    python run_methylation_pipeline.py -r ./raw_data -g ./genome.fa -o ./results \\
        --sample-groups-file ./sample_groups.txt -e
        
    # 完整参数示例 | Complete parameter example
    python run_methylation_pipeline.py \\
        -r ./raw_data -g ./genome.fa -o ./methylation_results \\
        -p ./target_promoter.fa -a ./annotation.gff3 \\
        -j 16 -c 8 -n 8 -d 0.25 -v 0.005 \\
        -w 200 -s 100 -m 300 \\
        --gene-bins 80 --flank-bins 25 --flank-size 3000 \\
        -e --verbose

样品分组文件格式示例 | Sample Groups File Format Examples:
    
    🔬 批量分组比较功能 | Batch Group Comparison Feature:
    - 支持多个分组的自动两两比较 | Support automatic pairwise comparison for multiple groups
    - 智能避免重复比较 (A vs B = B vs A) | Smart avoidance of duplicate comparisons  
    - 为每个比较生成独立结果文件 | Generate independent result files for each comparison
    
    方式1 (推荐) - sample_groups.txt (制表符分隔):
    # 基础两分组示例 | Basic two-group example
    FZY4201    CK
    FZY4202    Treatment  
    FZY4203    CK
    FZY4204    Treatment
    ➡️ 将进行: CK vs Treatment (1个比较)
    
    # 多分组示例 | Multi-group example
    FZY4201    Control
    FZY4202    Treatment_A
    FZY4203    Treatment_B
    FZY4204    Control
    FZY4205    Treatment_A
    FZY4206    Treatment_B
    ➡️ 将进行: Control vs Treatment_A, Control vs Treatment_B, Treatment_A vs Treatment_B (3个比较)
    
    # 复杂多分组示例 | Complex multi-group example  
    Sample1    Healthy
    Sample2    Disease_A
    Sample3    Disease_B
    Sample4    Disease_C
    Sample5    Healthy
    ➡️ 将进行: Healthy vs Disease_A, Healthy vs Disease_B, Healthy vs Disease_C, 
              Disease_A vs Disease_B, Disease_A vs Disease_C, Disease_B vs Disease_C (6个比较)
    
    方式2 - sample_groups.json:
    {
        "FZY4201": "Control",
        "FZY4202": "Treatment_A",
        "FZY4203": "Treatment_B", 
        "FZY4204": "Control"
    }

快速参数说明 | Quick Parameter Guide:
    -r, --raw-dir           📁 原始数据目录
    -g, --genome            🧬 基因组FASTA文件  
    -o, --output-dir        📂 输出目录
    -t, --target-fa         🎯 目标序列文件
    -p, --target-promoter-fa 🎯 目标启动子序列文件
    -a, --annotation-gff    📊 基因注释GFF文件
    -e, --enhanced-mode     🔬 启用增强分析模式
    -j, --threads           🔧 线程数
    -c, --min-coverage      📈 最小覆盖度
    -n, --min-cytosines     🧬 最小胞嘧啶数
    -d, --methylation-diff-threshold  📊 甲基化差异阈值
    -v, --pvalue-threshold  📈 p值阈值
    -w, --window-size       🪟 窗口大小
    -s, --step-size         👣 步长
    -m, --merge-distance    🔗 合并距离
    --sample-groups-file    🧬 样品分组文件 (TXT或JSON格式)
    --dry-run               🔍 试运行模式
    --verbose               📝 详细输出模式

样品分组文件示例 | Sample Groups File Example:
    创建sample_groups.txt文件，内容如下:
    FZY4201    CK
    FZY4202    Treatment
    FZY4203    CK  
    FZY4204    Treatment
    # 可以添加注释行
    
🔬 批量分组比较结果文件 | Batch Group Comparison Result Files:
    当使用多分组时，将为每个比较生成独立的结果文件:
    - {Group1}_vs_{Group2}_differential_sites.txt     # 所有差异位点
    - {Group1}_vs_{Group2}_hypermethylated_sites.txt  # 高甲基化位点  
    - {Group1}_vs_{Group2}_hypomethylated_sites.txt   # 低甲基化位点
    - {Group1}_vs_{Group2}_summary.txt                # 比较统计摘要
    - methylkit_comprehensive_summary.txt             # 所有比较的综合报告
    
    示例: 如果有Control, Treatment_A, Treatment_B三个分组，将生成:
    ├── Control_vs_Treatment_A_differential_sites.txt
    ├── Control_vs_Treatment_B_differential_sites.txt  
    ├── Treatment_A_vs_Treatment_B_differential_sites.txt
    └── methylkit_comprehensive_summary.txt
"""

from biopytools.methylation_pipeline.main import main

if __name__ == "__main__":
    main()
