#!/usr/bin/env python3
"""
ANNOVAR注释运行脚本 | ANNOVAR Annotation Runner Script
这是一个简化的入口脚本，用于运行基因注释分析 | Simple entry script for gene annotation analysis

用法 | Usage:
    python run_annovar.py --gff3 annotation.gff3 --genome genome.fa --vcf variants.vcf --build-ver OV
    
示例 | Examples:
    # 基本注释 | Basic annotation
    python run_annovar.py --gff3 gene.gff3 --genome genome.fa --vcf variants.vcf --build-ver OV
    
    # 指定数据库和输出路径 | Specify database and output paths
    python run_annovar.py --gff3 gene.gff3 --genome genome.fa --vcf variants.vcf --build-ver OV \\
        --annovar-path /path/to/annovar --database-path /path/to/db --output-dir results
    
    # 只运行特定步骤 | Run specific step only
    python run_annovar.py --gff3 gene.gff3 --genome genome.fa --vcf variants.vcf --build-ver OV --step 1
    
    # 启用VCF质量过滤 | Enable VCF quality filtering
    python run_annovar.py --gff3 gene.gff3 --genome genome.fa --vcf variants.vcf --build-ver OV \\
        --enable-vcf-filter --qual-threshold 30
"""

from biopytools.annovar.main import main

if __name__ == "__main__":
    main()
