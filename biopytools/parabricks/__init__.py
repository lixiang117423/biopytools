"""
🧬 parabricks WGS批处理分析工具包 | parabricks WGS Batch Processing Analysis Toolkit 🧬
功能 🚀: 全基因组测序数据的parabricks分析流程，支持批量处理和质量控制 | 
Features ✨: parabricks analysis pipeline for whole genome sequencing data with batch processing and quality control
作者 | Author: Xiang LI  
版本 🏷️ | Version: v1.0 - 模块化版本 | Modular version
日期 📅 | Date: 2025-08-08

使用示例 💡 | Usage Examples:
    from biopytools.parabricks import parabricksAnalyzer, parabricksConfig
    
    # 创建分析器 🛠️ | Create analyzer
    analyzer = parabricksAnalyzer(
        input_dir="/path/to/clean/data",
        output_dir="/path/to/output",
        reference="/path/to/reference.fa",
        threads=88
    )
    
    # 运行分析 🏃‍♂️ | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import parabricksAnalyzer
from .config import parabricksConfig

__all__ = ['parabricksAnalyzer', 'parabricksConfig']
