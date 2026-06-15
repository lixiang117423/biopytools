"""
HiFi+Hi-C基因组组装与挂载流程模块|HiFi+Hi-C Genome Assembly and Scaffolding Workflow Module

该模块整合了4个子模块，实现完整的植物基因组组装流程：
1. hifi_hic: HiFi组装（可选NGS polish）
2. haphic: Hi-C染色体挂载
3. rename_chromosomes: 染色体标准化命名
4. hic_heatmap: 全基因组Hi-C热图生成

This module integrates 4 sub-modules for complete plant genome assembly workflow:
1. hifi_hic: HiFi assembly (with optional NGS polish)
2. haphic: Hi-C chromosome scaffolding
3. rename_chromosomes: Chromosome standard naming
4. hic_heatmap: Whole genome Hi-C heatmap generation
"""

from .main import HifiHicWorkflow
from .config import HifiHicWorkflowConfig

__version__ = "1.0.0"
__all__ = ["HifiHicWorkflow", "HifiHicWorkflowConfig"]
