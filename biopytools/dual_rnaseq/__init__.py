"""
互作转录组分析工具包|Dual RNA-seq Analysis Toolkit
功能: 双物种转录组同时分析，包括物种分类、定量和表达矩阵生成|
Features: Simultaneous analysis of dual-species transcriptome, including species classification, quantification and expression matrix generation
应用场景|Applications:
  - 病原体-宿主互作|Pathogen-host interaction
  - 植物-微生物互作|Plant-microbe interaction
  - 共生体系研究|Symbiosis system research
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-09

使用示例|Usage Examples:
    from biopytools.dual_rnaseq import DualRNASeqAnalyzer, DualRNASeqConfig

    # 创建分析器|Create analyzer
    analyzer = DualRNASeqAnalyzer(
        species1_name="host",
        species1_genome="host.fa",
        species1_gtf="host.gtf",
        species2_name="pathogen",
        species2_genome="pathogen.fa",
        species2_gtf="pathogen.gtf",
        input_path="./fastq_data",
        output_dir="dual_rnaseq_results"
    )

    # 运行分析|Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import DualRNASeqAnalyzer
from .config import DualRNASeqConfig

__all__ = ['DualRNASeqAnalyzer', 'DualRNASeqConfig']
