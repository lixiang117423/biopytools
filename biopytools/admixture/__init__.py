"""
ADMIXTURE群体结构分析工具包 | ADMIXTURE Population Structure Analysis Toolkit
功能: VCF到ADMIXTURE分析的完整流程，支持群体结构分析和协变量生成 | 
Features: Complete pipeline from VCF to ADMIXTURE analysis, supporting population structure analysis and covariate generation
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化重构版 | Modular refactored version
日期 | Date: 2025-07-17

使用示例 | Usage Examples:
    from biopytools.admixture_toolkit import AdmixtureAnalyzer, AdmixtureConfig
    
    # 创建分析器 | Create analyzer
    analyzer = AdmixtureAnalyzer(
        vcf_file="data.vcf.gz",
        output_dir="admixture_results",
        min_k=2,
        max_k=10
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import AdmixtureAnalyzer
from .config import AdmixtureConfig

__all__ = ['AdmixtureAnalyzer', 'AdmixtureConfig']
