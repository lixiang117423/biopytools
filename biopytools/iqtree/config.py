"""
🔧 IQ-TREE 分析配置管理模块 | IQ-TREE Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

@dataclass
class TreeConfig:
    """🌲 系统发育树分析配置类 | Phylogenetic Tree Analysis Configuration Class"""
    
    # 输入文件 | Input files
    input_file: str
    output_dir: str = './iqtree_output'
    prefix: str = 'phylo_tree'
    
    # 核心参数 | Core parameters
    model: Optional[str] = None  # None表示自动选择 | None means auto-selection
    threads: int = 88
    
    # Bootstrap参数 | Bootstrap parameters
    bootstrap: int = 1000  # UFBoot重复次数 | UFBoot replicates
    bootstrap_type: str = 'ufboot'  # 'ufboot' or 'standard'
    save_boot_trees: bool = False  # 保存所有bootstrap树 | Save all bootstrap trees
    
    # 外群设置 | Outgroup settings
    outgroup: Optional[str] = None
    
    # 高级功能开关 | Advanced features switches
    partition_file: Optional[str] = None  # 分区文件 | Partition file
    partition_mode: str = 'edge-linked'  # 'edge-linked', 'edge-equal', 'edge-unlinked'
    enable_concordance: bool = False  # 一致性因子分析 | Concordance factor analysis
    concordance_trees: Optional[str] = None  # 用于gCF的树文件 | Tree file for gCF
    enable_ancestral: bool = False  # 祖先状态重建 | Ancestral state reconstruction
    constraint_tree: Optional[str] = None  # 约束树 | Constraint tree
    
    # 其他参数 | Other parameters
    seed: Optional[int] = None
    runs: int = 1
    redo: bool = False
    
    # 工具路径 | Tool path
    iqtree_path: str = 'iqtree'
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        if self.partition_file:
            self.partition_file = os.path.normpath(os.path.abspath(self.partition_file))
        
        if self.concordance_trees:
            self.concordance_trees = os.path.normpath(os.path.abspath(self.concordance_trees))
            
        if self.constraint_tree:
            self.constraint_tree = os.path.normpath(os.path.abspath(self.constraint_tree))
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查输入文件 | Check input file
        if not os.path.exists(self.input_file):
            errors.append(f"❌ 输入文件不存在 | Input file does not exist: {self.input_file}")
        
        # 检查分区文件 | Check partition file
        if self.partition_file and not os.path.exists(self.partition_file):
            errors.append(f"❌ 分区文件不存在 | Partition file does not exist: {self.partition_file}")
        
        # 检查一致性因子分析文件 | Check concordance factor files
        if self.enable_concordance:
            if not self.concordance_trees:
                errors.append(f"❌ 启用一致性因子分析需要提供树文件 | Concordance factor analysis requires tree file")
            elif not os.path.exists(self.concordance_trees):
                errors.append(f"❌ 一致性因子树文件不存在 | Concordance tree file does not exist: {self.concordance_trees}")
        
        # 检查约束树 | Check constraint tree
        if self.constraint_tree and not os.path.exists(self.constraint_tree):
            errors.append(f"❌ 约束树文件不存在 | Constraint tree file does not exist: {self.constraint_tree}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Thread number must be positive: {self.threads}")
        
        if self.bootstrap < 0:
            errors.append(f"❌ Bootstrap次数不能为负数 | Bootstrap replicates cannot be negative: {self.bootstrap}")
        
        if self.bootstrap_type not in ['ufboot', 'standard']:
            errors.append(f"❌ Bootstrap类型必须是'ufboot'或'standard' | Bootstrap type must be 'ufboot' or 'standard': {self.bootstrap_type}")
        
        if self.partition_mode not in ['edge-linked', 'edge-equal', 'edge-unlinked']:
            errors.append(f"❌ 分区模式无效 | Invalid partition mode: {self.partition_mode}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
