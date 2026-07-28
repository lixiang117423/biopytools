"""PSVCP 线性泛基因组构建|PSVCP linear pangenome construction

封装 PSVCP(psvcp_v1.01)的泛基因组构建流程:基于多个组装基因组,
用 MUMmer + Assemblytics 检测 PAV 并依次并入参考,构建线性泛基因组。
|Wraps PSVCP's pangenome construction: incorporate PAVs from each query
genome (MUMmer + Assemblytics) into the reference to build a linear pan-genome.

注意|Note: 不导出 main 函数,避免与子模块名冲突(参见 [[phobius-module]] 记忆)
|Do not export main() to avoid name collision with submodule name.
"""

from .config import PSVCPConfig
from .pipeline import PangenomeConstructor

__version__ = "1.0.0"

__all__ = [
    'PSVCPConfig',
    'PangenomeConstructor',
]
