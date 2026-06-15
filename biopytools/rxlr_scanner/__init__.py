"""
RxLR效应蛋白扫描工具包|RxLR Effector Protein Scanner Toolkit
功能: 批量扫描蛋白质序列中的RxLR和EER基序|
Features: Batch scanning for RxLR and EER motifs in protein sequences
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-05

使用示例|Usage Examples:
    from biopytools.rxlr_scanner import RxLRScanner, RxLRConfig

    # 创建扫描器|Create scanner
    scanner = RxLRScanner(
        input_file="proteins.fa",
        output_prefix="rxlr_results"
    )

    # 运行扫描|Run scanning
    scanner.run_scan()

    # 获取结果|Get results
    results = scanner.get_results()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import RxLRScanner
from .config import RxLRConfig
from .scanner import RxLRMotifScanner

__all__ = ['RxLRScanner', 'RxLRConfig', 'RxLRMotifScanner']
