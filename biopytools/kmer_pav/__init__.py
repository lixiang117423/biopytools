"""
🧬 K-mer PAV (Presence/Absence Variation) 分析工具包 | K-mer PAV Analysis Toolkit
功能: 从基因组提取k-mer并检查在fastq样本中的存在情况，生成存在性矩阵 | 
Features: Extract k-mers from genome and check their presence in fastq samples, generate presence matrix
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-08-05

使用示例 | Usage Examples:
    from biopytools.kmer_pav import KmerPAVAnalyzer, PAVConfig
    
    # 🚀 创建分析器 | Create analyzer
    analyzer = KmerPAVAnalyzer(
        genome_file="genome_chr12.fa",
        fastq_dir="/path/to/fastq",
        output_dir="kmer_pav_results",
        kmer_size=51
    )
    
    # 🏃 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import KmerPAVAnalyzer
from .config import PAVConfig

__all__ = ['KmerPAVAnalyzer', 'PAVConfig']
