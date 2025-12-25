"""
🧬 转录组预测分析工具包 | Transcriptome-based Prediction Analysis Toolkit
功能: RNA-seq数据的转录组预测完整流程，支持比对、组装、注释和编码区预测 | 
Features: Complete pipeline for transcriptome-based prediction from RNA-seq data, supporting alignment, assembly, annotation and CDS prediction
作者 | Author: Xiang LI  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-08-20

使用示例 | Usage Examples:
    from biopytools.transcriptome_prediction import TranscriptomeAnalyzer, TranscriptomeConfig
    
    # 创建分析器 | Create analyzer
    analyzer = TranscriptomeAnalyzer(
        genome_file="genome.fa",
        rna_seq_files=["sample1_R1.fq", "sample1_R2.fq"],
        output_dir="transcriptome_results"
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import TranscriptomeAnalyzer
from .config import TranscriptomeConfig

__all__ = ['TranscriptomeAnalyzer', 'TranscriptomeConfig']
