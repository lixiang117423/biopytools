"""
Cactus泛基因组分析模块|Cactus Pangenome Analysis Module

本模块提供了Cactus泛基因组构建流程的Python封装
This module provides Python wrapper for Cactus pangenome construction pipeline

主要功能|Main features:
    - 基于Minigraph-Cactus流程构建泛基因组|Build pangenome using Minigraph-Cactus pipeline
    - 支持Singularity容器调用|Support Singularity container invocation
    - 自动化流程管理|Automated pipeline management
    - 断点续传支持|Checkpoint resume support

作者|Author: Xiang LI
版本|Version: 1.0.0
"""

from .config import CactusConfig
from .pangenome import CactusPangenomeRunner

__version__ = "1.0.0"
__all__ = ["CactusConfig", "CactusPangenomeRunner"]
