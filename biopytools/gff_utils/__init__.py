"""
GFF3基因转录本提取工具包 | GFF3 Gene Transcript Extraction Toolkit
功能: 从GFF3文件中提取基因和转录本的整合信息 | 
Features: Extract integrated gene and transcript information from GFF3 files
作者 | Author: Xiang LI  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-07-10

使用示例 | Usage Examples:
    from biopytools.biohelpers.gff_utils import GFFAnalyzer, GFFConfig
    
    # 创建分析器 | Create analyzer
    analyzer = GFFAnalyzer(
        gff3_file="annotation.gff3",
        output_file="gene_transcript_info.tsv"
    )
    
    # 运行提取 | Run extraction
    analyzer.run_extraction()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import GFFAnalyzer
from .config import GFFConfig

__all__ = ['GFFAnalyzer', 'GFFConfig']
