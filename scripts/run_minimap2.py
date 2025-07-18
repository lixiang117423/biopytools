#!/usr/bin/env python3
"""
Minimap2比对运行脚本 | Minimap2 Alignment Runner Script
这是一个简化的入口脚本，用于运行Minimap2比对和未比对区间分析 | Simple entry script for Minimap2 alignment and unmapped region analysis

用法 | Usage:
    python run_minimap2.py -t target.fa -q query.fa -o minimap2_output --tp-type P
    
示例 | Examples:
    # 基本分析（默认使用tp:A:P类型）
    python run_minimap2.py -t reference.fa -q assembly.fa
    
    # 指定参数（短参数） | Specify parameters (short options)
    python run_minimap2.py -t ref.fa -q asm.fa -x asm10 -p 88 -m 2000 -u 5000 --tp-type P
    
    # 指定参数（长参数） | Specify parameters (long options)
    python run_minimap2.py --target ref.fa --query asm.fa --preset asm10 --threads 88 \\
        --min-match 2000 --min-unmapped 5000 --tp-type SP
    
    # 使用不同预设和tp类型 | Use different presets and tp types
    python run_minimap2.py -t ref.fa -q long_reads.fa -x map-ont -p 16 --tp-type S
    
    # 完整示例（对应您的原始命令） | Complete example (corresponding to your original command)
    python run_minimap2.py \\
        -t /share/org/YZWL/yzwl_lixg/project/08.rice_pangenome/01.data/all.genome.fa \\
        -q /share/org/YZWL/yzwl_lixg/project/08.rice_pangenome/01.data/all.genome.fa \\
        -o /share/org/YZWL/yzwl_lixg/project/08.rice_pangenome/02.minimap \\
        -x asm5 -p 88 --tp-type P
"""

from biopytools.minimap2.main import main

if __name__ == "__main__":
    main()
