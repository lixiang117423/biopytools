"""
多序列比对工具包 | Multiple Sequence Alignment Toolkit
功能: 支持多种MSA软件的统一接口，包含自动化分析流程 | 
Features: Unified interface for multiple MSA tools with automated analysis pipeline
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-10-07

使用示例 | Usage Examples:
    from biopytools.msa_toolkit import MSAAnalyzer, MSAConfig
    
    # 创建分析器 | Create analyzer
    analyzer = MSAAnalyzer(
        input_file="sequences.fasta",
        output_prefix="alignment",
        method="mafft",
        threads=88
    )
    
    # 运行比对 | Run alignment
    analyzer.run_alignment()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import MSAAnalyzer
from .config import MSAConfig

__all__ = ['MSAAnalyzer', 'MSAConfig']
