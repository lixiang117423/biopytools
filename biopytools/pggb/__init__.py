"""
PGGB泛基因组图构建模块|PGGB Pangenome Graph Builder Module

本模块提供了PGGB泛基因组图构建流程的Python封装
This module provides Python wrapper for PGGB pangenome graph construction pipeline

主要功能|Main features:
    - 基于wfmash+seqwish+smoothxg构建泛基因组变异图|Build pangenome variation graph using wfmash+seqwish+smoothxg
    - 支持conda环境自动检测与调用|Support conda environment auto-detection and invocation
    - 自动生成fasta索引|Auto-generate fasta index
    - 断点续传支持(通过pggb -r)|Resume support (via pggb -r)

作者|Author: Xiang LI
版本|Version: 1.0.0
"""

from .main import PGGBRunner
from .config import PGGBConfig

__version__ = "1.0.0"
__all__ = ["PGGBRunner", "PGGBConfig"]
