"""
多序列比对配置管理模块 | MSA Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class MSAConfig:
    """MSA配置类 | MSA Configuration Class"""
    
    # 输入输出 | Input/Output
    input_file: str
    output_prefix: str = 'alignment'
    
    # 比对方法 | Alignment method
    method: str = 'mafft'  # mafft, clustalo, muscle, t_coffee
    
    # 通用参数 | Common parameters
    threads: int = 88
    
    # MAFFT特定参数 | MAFFT-specific parameters
    mafft_strategy: str = 'auto'  # auto, linsi, ginsi, einsi, fftns, fftnsi
    mafft_maxiterate: int = 1000
    
    # Clustal Omega特定参数 | Clustal Omega-specific parameters
    clustalo_iterations: int = 0  # 0表示不迭代 | 0 means no iteration
    
    # MUSCLE特定参数 | MUSCLE-specific parameters  
    muscle_maxiters: int = 16
    
    # 输出格式 | Output format
    output_format: str = 'fasta'  # fasta, clustal, phylip, nexus
    
    # 工具路径 | Tool paths
    mafft_path: str = 'mafft'
    clustalo_path: str = 'clustalo'
    muscle_path: str = 'muscle'
    tcoffee_path: str = 't_coffee'
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 标准化路径 | Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        
        # 验证方法 | Validate method
        valid_methods = ['mafft', 'clustalo', 'muscle', 't_coffee']
        if self.method.lower() not in valid_methods:
            raise ValueError(f"❌ 不支持的比对方法 | Unsupported method: {self.method}. "
                           f"可选 | Available: {', '.join(valid_methods)}")
        
        self.method = self.method.lower()
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查输入文件 | Check input file
        if not os.path.exists(self.input_file):
            errors.append(f"❌ 输入文件不存在 | Input file does not exist: {self.input_file}")
        
        # 检查线程数 | Check thread count
        if self.threads < 1:
            errors.append(f"❌ 线程数必须≥1 | Thread count must be ≥1: {self.threads}")
        
        # 检查MAFFT策略 | Check MAFFT strategy
        if self.method == 'mafft':
            valid_strategies = ['auto', 'linsi', 'ginsi', 'einsi', 'fftns', 'fftnsi']
            if self.mafft_strategy not in valid_strategies:
                errors.append(f"❌ 无效的MAFFT策略 | Invalid MAFFT strategy: {self.mafft_strategy}. "
                            f"可选 | Available: {', '.join(valid_strategies)}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
