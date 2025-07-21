#!/usr/bin/env python3
"""
ADMIXTURE群体结构分析运行脚本 | ADMIXTURE Population Structure Analysis Runner Script
这是一个简化的入口脚本，用于运行ADMIXTURE分析 | Simple entry script for running ADMIXTURE analysis

用法 | Usage:
    # 使用长参数
    python run_admixture.py --vcf data.vcf.gz --output admixture_results
    
    # 使用短参数
    python run_admixture.py -v data.vcf.gz -o admixture_results -t 8
    
示例 | Examples:
    # 基本分析 | Basic analysis
    python run_admixture.py -v final_filtered.vcf.gz -o results
    
    # 指定K值范围和线程数 | Specify K range and thread count
    python run_admixture.py -v data.vcf.gz -o results -k 2 -K 8 -t 8
    
    # 跳过预处理步骤 | Skip preprocessing steps
    python run_admixture.py -v clean_data.vcf.gz -o results -s
    
    # 自定义质控参数 | Custom QC parameters
    python run_admixture.py -v data.vcf.gz -o results -m 0.05 -M 0.05 -H 1e-5
    
    # 完整参数示例 | Complete parameter example
    python run_admixture.py -v data.vcf.gz -o results \\
        -k 2 -K 10 -c 5 -t 8 -m 0.01 -M 0.1 -H 1e-6 -i
    
参数说明 | Parameter Description:
    -v, --vcf        输入VCF文件 | Input VCF file
    -o, --output     输出目录 | Output directory
    -k, --min-k      最小K值 | Minimum K value
    -K, --max-k      最大K值 | Maximum K value
    -c, --cv-folds   交叉验证折数 | Cross-validation folds
    -t, --threads    线程数 | Number of threads
    -m, --maf        MAF阈值 | MAF threshold
    -M, --missing    缺失率阈值 | Missing rate threshold
    -H, --hwe        HWE p值阈值 | HWE p-value threshold
    -s, --skip-preprocessing  跳过预处理 | Skip preprocessing
    -i, --keep-intermediate   保留中间文件 | Keep intermediate files

快速使用指南 | Quick Usage Guide:
    # 最简单的使用方法 | Simplest usage
    python run_admixture.py -v my_data.vcf.gz -o results
    
    # 快速分析（较小K值范围） | Quick analysis (smaller K range)
    python run_admixture.py -v my_data.vcf.gz -o results -k 2 -K 5 -t 8
    
    # 完整分析（更严格的质控） | Full analysis (stricter QC)
    python run_admixture.py -v my_data.vcf.gz -o results -k 2 -K 10 \\
        -m 0.05 -M 0.05 -H 1e-5 -t 16 -c 10
    
    # 高质量数据（跳过预处理） | High quality data (skip preprocessing)
    python run_admixture.py -v clean_data.vcf.gz -o results -s -t 16
    
    # 调试模式（保留中间文件） | Debug mode (keep intermediate files)
    python run_admixture.py -v my_data.vcf.gz -o results -i
"""

from biopytools.admixture.main import main

if __name__ == "__main__":
    main()
