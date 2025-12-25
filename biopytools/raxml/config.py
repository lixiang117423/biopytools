"""
🌳 RAxML分析配置管理模块 | RAxML Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List

@dataclass
class RAxMLConfig:
    """RAxML分析配置类 | RAxML Analysis Configuration Class"""
    
    # 必需参数 | Required parameters
    sequence_file: str
    output_name: str
    model: str = "GTRGAMMA"
    
    # 输出设置 | Output settings
    output_dir: str = './raxml_output'
    
    # 算法参数 | Algorithm parameters
    algorithm: str = "d"  # -f parameter (default: rapid hill-climbing)
    parsimony_seed: Optional[int] = None  # -p parameter
    bootstrap_seed: Optional[int] = None  # -b parameter
    rapid_bootstrap_seed: Optional[int] = None  # -x parameter
    runs: str = "1"  # -# parameter
    
    # 树搜索参数 | Tree search parameters
    starting_tree: Optional[str] = None  # -t parameter
    constraint_tree: Optional[str] = None  # -g parameter
    outgroup: Optional[str] = None  # -o parameter
    
    # 模型参数 | Model parameters
    categories: int = 25  # -c parameter
    rate_het_model: bool = True  # -V parameter (disable with False)
    gamma_median: bool = False  # -u parameter
    invariant_sites: bool = False  # Include invariant sites estimate
    
    # Bootstrap参数 | Bootstrap parameters
    bootstrap_convergence: Optional[str] = None  # -I parameter (autoFC|autoMR|autoMRE|autoMRE_IGN)
    bootstop_threshold: float = 0.03  # -B parameter
    bootstop_perms: int = 100  # --bootstop-perms parameter
    print_bootstrap_trees: bool = False  # -k parameter
    
    # 性能参数 | Performance parameters
    threads: int = 88  # -T parameter (PTHREADS version only)
    memory_saving: bool = False  # -U parameter
    pattern_compression: bool = True  # -H parameter (disable with False)
    
    # 高级参数 | Advanced parameters
    likelihood_epsilon: float = 0.1  # -e parameter
    ml_search_convergence: bool = False  # -D parameter
    random_starting_tree: bool = False  # -d parameter
    
    # 工具路径 | Tool paths
    raxml_path: str = 'raxmlHPC-PTHREADS'
    
    # 质量控制 | Quality control
    no_seq_check: bool = False  # --no-seq-check
    silent: bool = False  # --silent
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.sequence_file = os.path.normpath(os.path.abspath(self.sequence_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        if self.starting_tree:
            self.starting_tree = os.path.normpath(os.path.abspath(self.starting_tree))
        
        if self.constraint_tree:
            self.constraint_tree = os.path.normpath(os.path.abspath(self.constraint_tree))
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查序列文件 | Check sequence file
        if not os.path.exists(self.sequence_file):
            errors.append(f"序列文件不存在 | Sequence file does not exist: {self.sequence_file}")
        
        # 检查起始树文件 | Check starting tree file
        if self.starting_tree and not os.path.exists(self.starting_tree):
            errors.append(f"起始树文件不存在 | Starting tree file does not exist: {self.starting_tree}")
        
        # 检查约束树文件 | Check constraint tree file
        if self.constraint_tree and not os.path.exists(self.constraint_tree):
            errors.append(f"约束树文件不存在 | Constraint tree file does not exist: {self.constraint_tree}")
        
        # 检查参数范围 | Check parameter ranges
        if self.categories <= 0:
            errors.append(f"分类数量必须为正整数 | Number of categories must be positive: {self.categories}")
        
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Number of threads must be positive: {self.threads}")
        
        if not 0 < self.bootstop_threshold < 1:
            errors.append(f"Bootstrap停止阈值必须在0-1之间 | Bootstrap stop threshold must be between 0-1: {self.bootstop_threshold}")
        
        # 检查模型格式 | Check model format
        valid_models = [
            "GTRCAT", "GTRCATI", "GTRGAMMA", "GTRGAMMAI",
            "PROTCAT", "PROTGAMMA", "PROTGAMMAI",
            "BINCAT", "BINGAMMA", "BINGAMMAI",
            "MULTICAT", "MULTIGAMMA", "MULTIGAMMAI"
        ]
        
        model_base = self.model.replace("X", "").replace("F", "")
        if not any(model_base.startswith(vm) for vm in valid_models):
            errors.append(f"不支持的模型 | Unsupported model: {self.model}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
