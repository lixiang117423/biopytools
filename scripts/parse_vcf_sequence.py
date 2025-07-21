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

# import sys

# def main_wrapper():
#     """主函数包装器，确保正确导入和调用"""
#     try:
#         # 尝试导入main函数
#         from biopytools.vcf_sequence_toolkit.main import main
#         if callable(main):
#             return main()
#         else:
#             raise ImportError("main is not callable")
#     except ImportError as e:
#         try:
#             # 备用方案：直接导入模块并调用
#             from biopytools.vcf_sequence_toolkit import main as main_module
#             if hasattr(main_module, 'main') and callable(main_module.main):
#                 return main_module.main()
#             else:
#                 raise ImportError("Cannot find callable main function")
#         except ImportError:
#             print(f"Error importing sequence_toolkit: {e}")
#             print("请确保已正确安装vcf_sequence_toolkit模块 | Please ensure vcf_sequence_toolkit module is properly installed")
#             sys.exit(1)

# if __name__ == "__main__":
#     main_wrapper()

from biopytools.vcf_sequence_toolkit.main import main

if __name__ == "__main__":
    main()
