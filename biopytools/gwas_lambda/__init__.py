"""
GWAS Lambda计算工具 | GWAS Lambda GC Calculator
功能: 批量分析GWAS结果，计算Lambda GC值并评估群体分层情况 |
Features: Batch analysis of GWAS results, calculate Lambda GC values and assess population stratification
作者 | Author: Xiang LI
版本 | Version: 1.0.0 - 模块化版本 | Modular version
日期 | Date: 2025-12-19

使用示例 | Usage Examples:
    from biopytools.gwas_lambda import GWASLambdaCalculator, GWASLambdaConfig

    # 创建计算器 | Create calculator
    calculator = GWASLambdaCalculator(
        search_pattern="feture_*/GWAS_Result.mlm.manht_input",
        output_file="lambda_assessment.txt",
        significance_threshold=1e-5
    )

    # 运行分析 | Run analysis
    calculator.run_batch_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import GWASLambdaCalculator
from .config import GWASLambdaConfig
from .analyzer import GWASResultAnalyzer

__all__ = ['GWASLambdaCalculator', 'GWASLambdaConfig', 'GWASResultAnalyzer']