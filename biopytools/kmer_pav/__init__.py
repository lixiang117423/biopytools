"""
K-mer PAV (Presence/Absence Variation) 分析工具包 (双阶段设计)
功能: 基于KMC工具的双阶段k-mer变异检测，包含k-mer来源追踪
  阶段1: 从多个文件构建k-mer数据库
  阶段2: 样本查询和比较分析
  特性: 追踪k-mer在database和query中的分布
作者 | Author: Xiang LI  
版本 | Version: v3.1 - 双阶段设计版（添加from列）
日期 | Date: 2025-07-18

设计理念 | Design Philosophy:
  1. 数据库构建: 从大量文件中提取所有k-mer
  2. 样本分析: 查询文件与数据库的比较
  3. 灵活处理: FASTQ整文件，FASTA单序列
  4. 来源追踪: 显示k-mer来自database中的哪个文件

使用示例 | Usage Examples:
    from biopytools.kmer_pav import KmerPAVAnalyzer
    
    # 创建分析器 | Create analyzer
    analyzer = KmerPAVAnalyzer(
        database_input="database_files/",
        query_input="query_samples/",
        kmer_size=31,
        output_prefix="kmer_analysis"
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "3.1.0"
__author__ = "Xiang LI"

from .main import KmerPAVAnalyzer
from .config import KmerConfig

__all__ = ['KmerPAVAnalyzer', 'KmerConfig']
