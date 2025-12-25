"""
🧬 GTX WGS批处理分析工具包 | GTX WGS Batch Processing Analysis Toolkit 🧬
功能 🚀: 全基因组测序数据的GTX分析流程，支持批量处理和质量控制 | 
Features ✨: GTX analysis pipeline for whole genome sequencing data with batch processing and quality control
作者 | Author: Xiang LI  
版本 🏷️ | Version: v1.0 - 模块化版本 | Modular version
日期 📅 | Date: 2025-08-08

使用示例 💡 | Usage Examples:
    from biopytools.gtx import GTXAnalyzer, GTXConfig
    
    # 创建分析器 🛠️ | Create analyzer
    analyzer = GTXAnalyzer(
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

from .main import GTXAnalyzer
from .config import GTXConfig

__all__ = ['GTXAnalyzer', 'GTXConfig']
