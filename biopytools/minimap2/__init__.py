"""
Minimap2比对和未比对区间提取工具包 | Minimap2 Alignment and Unmapped Region Extraction Toolkit
功能: Minimap2全基因组比对、PAF结果解析和未比对区间提取的完整流程 | 
Features: Complete pipeline for Minimap2 whole genome alignment, PAF result parsing and unmapped region extraction
作者 | Author: Xiang LI  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-07-18

使用示例 | Usage Examples:
    from biopytools.minimap2 import Minimap2Analyzer, Minimap2Config
    
    # 创建分析器 | Create analyzer
    analyzer = Minimap2Analyzer(
        target_genome="target.fa",
        query_genome="query.fa",
        output_dir="minimap2_output",
        preset="asm5",
        tp_type="P"  # 可选：'S', 'P', 'SP'
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import Minimap2Analyzer
from .config import Minimap2Config

__all__ = ['Minimap2Analyzer', 'Minimap2Config']
