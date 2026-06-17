"""
Hi-C数据质量控制评估工具包|Hi-C Data Quality Control Assessment Toolkit
功能: 使用pairtools评估Hi-C mapping数据质量|
Features: Assess Hi-C mapping data quality using pairtools
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-10

使用示例|Usage Examples:
    from biopytools.pairtools_qc import PairtoolsQCRunner, PairtoolsQCConfig

    # 创建QC评估器|Create QC assessor
    qc = PairtoolsQCRunner(
        pairs_file="sample.pairs.gz",
        pairtools_path="~/miniforge3/envs/pairtools_v.1.1.3/bin/pairtools"
    )

    # 运行评估|Run assessment
    results = qc.run()
    qc.print_report(results)
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import PairtoolsQCRunner
from .config import PairtoolsQCConfig
from .calculator import PairtoolsQCCalculator

__all__ = ['PairtoolsQCRunner', 'PairtoolsQCConfig', 'PairtoolsQCCalculator']
