"""
RNA-seq分析工具包 | RNA-seq Analysis Toolkit
功能: HISAT2比对、StringTie定量和表达矩阵合并的完整流程 | 
Features: Complete pipeline for HISAT2 alignment, StringTie quantification and expression matrix merging
作者 | Author: Xiang LI  
版本 | Version: v10 - 模块化重构版 | Modular refactored version
日期 | Date: 2025-07-10

使用示例 | Usage Examples:
    from biopytools.rnaseq import RNASeqAnalyzer, RNASeqConfig
    
    # 创建分析器 | Create analyzer
    analyzer = RNASeqAnalyzer(
        genome_file="genome.fa",
        gtf_file="annotation.gtf",
        input_path="/path/to/fastq",
        output_dir="rnaseq_results"
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "10.0.0"
__author__ = "Xiang LI"

from .main import RNASeqAnalyzer
from .config import RNASeqConfig

__all__ = ['RNASeqAnalyzer', 'RNASeqConfig']
