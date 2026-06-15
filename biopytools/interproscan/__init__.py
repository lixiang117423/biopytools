"""
InterProScan蛋白质功能注释工具包|InterProScan Protein Function Annotation Toolkit
功能: 对蛋白质序列进行功能结构域注释和GO术语预测|
Features: Functional domain annotation and GO term prediction for protein sequences
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-01-12

使用示例|Usage Examples:
    from biopytools.interproscan import InterProScanAnnotator, InterProScanConfig

    # 创建注释器|Create annotator
    annotator = InterProScanAnnotator(
        input_file="proteins.fa",
        output_prefix="results",
        interproscan_path="/path/to/interproscan.sh"
    )

    # 运行注释|Run annotation
    annotator.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import InterProScanAnnotator
from .config import InterProScanConfig
from .parser import InterProScanParser, ProteinMatch, ProteinSummary
from .formatter import InterProScanFormatter
from .go_database import GODatabase

__all__ = [
    'InterProScanAnnotator',
    'InterProScanConfig',
    'InterProScanParser',
    'ProteinMatch',
    'ProteinSummary',
    'InterProScanFormatter',
    'GODatabase'
]
