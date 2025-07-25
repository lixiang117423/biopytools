#!/usr/bin/env python3
"""
基因组组装分析运行脚本 | Genome Assembly Analysis Runner Script
这是一个简化的入口脚本，用于运行基因组组装分析 | Simple entry script for running genome assembly analysis

用法 | Usage:
    python run_genome_assembly.py --hifi-reads hifi.fastq.gz -o assembly_results
    
示例 | Examples:
    # 基本组装分析 | Basic assembly analysis
    python run_genome_assembly.py --hifi-reads hifi.fastq.gz -o assembly_results
    
    # 包含ONT数据 | Include ONT data
    python run_genome_assembly.py --hifi-reads hifi.fq.gz --ont-reads ont.fq.gz --genome-size 3.2g
    
    # 家系模式分析 | Family trio analysis
    python run_genome_assembly.py --hifi-reads child_hifi.fq.gz --trio-mode \\
        --parent1-reads p1.fq.gz --parent2-reads p2.fq.gz --child-reads child.fq.gz
    
    # 高性能配置 | High-performance configuration
    python run_genome_assembly.py --hifi-reads hifi.fastq.gz \\
        --reference-genome T2T-CHM13.fa --threads 64 --memory 256
    
    # 自定义质控参数 | Custom quality control parameters
    python run_genome_assembly.py --hifi-reads hifi.fastq.gz \\
        --skip-contamination-screen --skip-flagger
    
    # 指定工具路径 | Specify tool paths
    python run_genome_assembly.py --hifi-reads hifi.fastq.gz \\
        --verkko-path /path/to/verkko --hifiasm-path /path/to/hifiasm

依赖工具 | Required Tools:
- Verkko (v.1.4.1+): 主要组装器 | Primary assembler
- hifiasm (v.0.19.6+): 超长读长组装器 | Ultra-long read assembler  
- Graphasing (v.0.3.1-alpha+): 分期管道 | Phasing pipeline
- NCBI FCS (v.0.4.0+): 外源污染筛查 | Foreign contamination screening
- Flagger (v.0.3.3+): 组装错误注释 | Assembly error annotation
- Merqury (v.1.0+): 质量评估 | Quality assessment
- NucFreq: 核苷酸频率分析 | Nucleotide frequency analysis
- Inspector (v.1.2+): 组装检查 | Assembly inspection
- DeepVariant (v.1.6+): 变异检测质量评估 | Variant calling quality assessment
- compleasm (v.0.2.5+): 基因完整性评估 | Gene completeness assessment
- minimap2 (v.2.26+): 序列比对 | Sequence alignment
- mashmap (v.3.1.3+): 快速序列比对 | Fast sequence mapping

安装方法 | Installation:
  参考各工具官方文档进行安装 | Refer to official documentation for each tool
  建议使用conda环境管理依赖 | Recommend using conda for dependency management
"""

from biopytools.genome_assembly.main import main

if __name__ == "__main__":
    main()