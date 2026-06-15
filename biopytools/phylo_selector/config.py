"""
系统发育树样品选择配置管理模块|Phylogenetic Tree Sample Selector Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class PhyloSelectorConfig:
    """系统发育树样品选择配置类|Phylogenetic Tree Sample Selector Configuration Class

    基于均匀间隔和PCA去重的智能选择算法
    Intelligent selection algorithm based on uniform interval and PCA deduplication
    """

    # 必需文件|Required files
    hierarchy_file: str  # 层级关系文件（必需，仅用于获取样品列表）|Hierarchy file (required, for sample list only)
    pca_file: str  # PCA坐标文件（必需）|PCA coordinates file (required)
    output_prefix: str  # 输出文件前缀|Output file prefix

    # 可选文件|Optional files
    newick_file: Optional[str] = None  # Newick树文件（可选，已弃用）|Newick tree file (optional, deprecated)
    group_file: Optional[str] = None  # 分组文件（可选，暂未使用）|Group file (optional, unused for now)

    # 选择参数|Selection parameters
    n_samples: int = 150  # 选择样品总数|Total number of samples to select
    min_samples_per_group: int = 1  # 每组最小样品数（暂未使用）|Minimum samples per group (unused for now)
    hierarchy_level: int = 10  # 层级深度（已弃用，保留以兼容）|Hierarchy level (deprecated, kept for compatibility)
    pca_dedup_threshold: float = 0.0001  # PCA去重阈值|PCA deduplication threshold

    # 输出控制|Output control
    generate_report: bool = True
    generate_csv: bool = True
    generate_visualization: bool = True

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.hierarchy_file = os.path.normpath(os.path.abspath(self.hierarchy_file))
        self.pca_file = os.path.normpath(os.path.abspath(self.pca_file))
        self.output_prefix = os.path.normpath(os.path.abspath(self.output_prefix))

        if self.newick_file:
            self.newick_file = os.path.normpath(os.path.abspath(self.newick_file))

        if self.group_file:
            self.group_file = os.path.normpath(os.path.abspath(self.group_file))

        # 创建输出目录|Create output directory
        self.output_dir = os.path.dirname(self.output_prefix)
        if self.output_dir:
            Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.hierarchy_file):
            errors.append(f"层级文件不存在|Hierarchy file does not exist: {self.hierarchy_file}")

        if not os.path.exists(self.pca_file):
            errors.append(f"PCA文件不存在|PCA file does not exist: {self.pca_file}")

        if self.newick_file and not os.path.exists(self.newick_file):
            errors.append(f"Newick文件不存在|Newick file does not exist: {self.newick_file}")

        if self.group_file and not os.path.exists(self.group_file):
            errors.append(f"分组文件不存在|Group file does not exist: {self.group_file}")

        # 检查参数范围|Check parameter ranges
        if self.n_samples <= 0:
            errors.append(f"样品数必须为正数|Number of samples must be positive: {self.n_samples}")

        if self.min_samples_per_group < 0:
            errors.append(f"每组最小样品数必须为非负数|Minimum samples per group must be non-negative: {self.min_samples_per_group}")

        if self.hierarchy_level <= 0:
            errors.append(f"层级深度必须为正数|Hierarchy level must be positive: {self.hierarchy_level}")

        if self.pca_dedup_threshold < 0:
            errors.append(f"PCA去重阈值必须为非负数|PCA dedup threshold must be non-negative: {self.pca_dedup_threshold}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
