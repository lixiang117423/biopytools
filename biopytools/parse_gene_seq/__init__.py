"""
基因序列提取工具包 🧬 | Gene Sequence Extraction Toolkit
功能: 从基因组FASTA文件和GFF注释文件中提取每个基因的DNA序列 |
Features: Extract gene DNA sequences from genome FASTA and GFF annotation files
作者 | Author: Gene Extractor Bot
版本 | Version: v1.0.0 - 模块化版本 | Modular version
日期 | Date: 2025-09-29

使用示例 | Usage Examples:
    from biopytools.parse_gene_seq import GeneSequenceExtractor, ExtractionConfig
    
    # 创建提取器 | Create extractor
    extractor = GeneSequenceExtractor(
        genome_file="genome.fasta",
        gff_file="annotation.gff",
        output_file="genes.fasta",
        feature_type="gene",
        threads=88
    )
    
    # 运行提取 | Run extraction
    extractor.run_extraction()
"""

__version__ = "1.0.0"
__author__ = "Gene Extractor Bot"

from .main import GeneSequenceExtractor
from .config import ExtractionConfig

__all__ = ['GeneSequenceExtractor', 'ExtractionConfig']
