"""
🧬 覆盖度分析工具包 | Depth Analysis Toolkit
功能: BAM/SAM文件覆盖度分析的完整流程，支持染色体和区间筛选、多线程并行处理 | 
Features: Complete pipeline for BAM/SAM depth analysis, supporting chromosome and interval filtering, multi-threading
作者 | Author: Xiang LI  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-08-22

使用示例 | Usage Examples:
    from biopytools.depth_analyzer import DepthAnalyzer, DepthConfig
    
    # 创建分析器 | Create analyzer
    analyzer = DepthAnalyzer(
        input_files=["sample1.bam", "sample2.bam"],
        output_file="depth_results.txt",
        chromosome="chr12",
        region="136491092:138554123",
        threads=88
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import DepthAnalyzer
from .config import DepthConfig

__all__ = ['DepthAnalyzer', 'DepthConfig']
