"""
通用K-mer分析工具包 | Universal K-mer Analysis Toolkit
支持FASTA/FASTQ文件的灵活k-mer分析，自动识别文件角色，高性能批量处理
Supports flexible k-mer analysis for FASTA/FASTQ files with automatic role detection and high-performance batch processing

主要功能 | Main Features:
- 智能文件格式识别和角色分配 | Intelligent file format detection and role assignment
- 基于KMC3的高性能k-mer计数 | High-performance k-mer counting based on KMC3
- 位置信息追踪和样品管理 | Position tracking and sample management
- 批量文件处理 | Batch file processing
- 多种输出格式支持 | Multiple output format support
- 滑窗分析 | Sliding window analysis

使用示例 | Usage Examples:
    from biopytools.kmer_universal import KmerAnalyzer, KmerConfig
    
    # 自动模式分析
    analyzer = KmerAnalyzer()
    analyzer.auto_analyze(
        input_paths=["/data/genomes/*.fa", "/data/samples/*.fq.gz"],
        output_dir="results/",
        kmer_size=51
    )
    
    # 明确指定角色
    analyzer.explicit_analyze(
        kmer_sources=["/ref/genome.fa"],
        query_targets=["/samples/*.fastq.gz"],
        output_dir="results/"
    )
"""

__version__ = "1.0.0"
__author__ = "BioPyTools Team"

from .main import KmerAnalyzer
from .config import KmerConfig
from .file_manager import FileManager, FileInfo
from .kmc_interface import KMCInterface

__all__ = ['KmerAnalyzer', 'KmerConfig', 'FileManager', 'FileInfo', 'KMCInterface']
