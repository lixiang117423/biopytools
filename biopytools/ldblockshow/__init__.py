"""
连锁不平衡热图分析工具包|LD Heatmap Analysis Toolkit
功能: 基于VCF文件生成连锁不平衡(LD)热图，支持GWAS和基因注释可视化|
Features: Generate linkage disequilibrium (LD) heatmap from VCF files, support GWAS and gene annotation visualization
应用场景|Applications:
  - 连锁不平衡分析|Linkage disequilibrium analysis
  - GWAS结果可视化|GWAS result visualization
  - 基因组区域LD结构分析|Genomic region LD structure analysis
作者|Author: hewm2008 (Original), Xiang LI (Python wrapper)
版本|Version: 1.0.0
日期|Date: 2026-02-09

使用示例|Usage Examples:
    from biopytools.ldblockshow import LDBlockShowAnalyzer, LDBlockShowConfig

    # 创建分析器|Create analyzer
    analyzer = LDBlockShowAnalyzer(
        vcf_file="variants.vcf.gz",
        output_prefix="ld_result",
        region="chr1:1000000-2000000"
    )

    # 运行分析|Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import LDBlockShowAnalyzer
from .config import LDBlockShowConfig

__all__ = ['LDBlockShowAnalyzer', 'LDBlockShowConfig']
