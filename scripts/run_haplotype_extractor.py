#!/usr/bin/env python3
"""
VCF单体型提取运行脚本 v2.0 | VCF Haplotype Extractor Runner Script v2.0
这是一个简化的入口脚本，用于运行VCF单体型提取分析 | Simple entry script for running VCF haplotype extraction analysis

用法 | Usage:
    python run_haplotype_extractor.py -v variants.vcf -p positions.txt -o haplotypes.txt
    
示例 | Examples:
    # 使用位置文件批量提取 | Batch extraction using position file
    python run_haplotype_extractor.py -v input.vcf.gz -p positions.txt -o output.txt
    
    # 提取单个SNP | Extract single SNP
    python run_haplotype_extractor.py -v input.vcf.gz --single -c OV12 -s 93635286 -o output.txt
    
    # 指定bcftools路径 | Specify bcftools path
    python run_haplotype_extractor.py -v input.vcf.gz -p positions.txt -o output.txt --bcftools-path /path/to/bcftools

位置文件格式 | Position File Format:
    # 可以有表头（推荐）| Can have header (recommended)
    CHR	POS
    OV12	93635286
    OV12	93440535
    OV12	93234780
    
    # 也可以没有表头 | Can also without header
    OV12	93635286
    OV12	93440535
    OV12	93234780

输出格式 | Output Format:
    CHROM	POS	REF	ALT	Sample1	Sample2	Sample3
    OV12	93635286	C	A	CA	AA	AA
    OV12	93440535	C	A	CA	AA	CA

更新日志 | Changelog:
    v2.0: 支持位置文件批量提取和标准矩阵输出格式
    v1.0: 初始版本，支持单个位点提取
"""

from biopytools.haplotype_extractor.main import main

if __name__ == "__main__":
    main()
