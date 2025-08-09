"""
K-mer分析配置管理模块 | K-mer Analysis Configuration Management Module ⚙️🛠️
"""

import os
import yaml
from pathlib import Path
from typing import List, Optional, Dict, Any
from dataclasses import dataclass, field
from enum import Enum

class FileRole(Enum):
    """文件角色枚举 🏷️"""
    KMER_SOURCE = "kmer_source"      # k-mer库来源 📚
    QUERY_TARGET = "query_target"    # 查询目标 🎯
    AUTO_DETECT = "auto_detect"      # 自动检测 ✨

class AssignmentStrategy(Enum):
    """角色分配策略枚举 🤔"""
    EXPLICIT = "explicit"            # 明确指定 👉
    SIZE_BASED = "size_based"        # 基于文件大小 ⚖️
    TYPE_BASED = "type_based"        # 基于文件类型 🏷️
    INTELLIGENT = "intelligent"      # 智能混合策略 🧠
    INTERACTIVE = "interactive"      # 交互式选择 💬

@dataclass
class KmerConfig:
    """K-mer分析配置类 ⚙️"""
    
    # 基本参数 | Basic Parameters 🔧
    kmer_size: int = 51
    threads: int = 88
    memory_gb: int = 880
    
    # KMC参数 | KMC Parameters 🔢
    min_count: int = 1              # -ci参数，最小计数 📉
    max_count: int = int(1e9)       # -cx参数，最大计数 📈
    canonical_form: bool = True     # 是否转换为标准形式 🔄
    signature_length: int = 9       # -p参数，签名长度 ✍️
    ram_only_mode: bool = False     # -r参数，仅RAM模式 ⚡
    strict_memory: bool = True      # -sm参数，严格内存模式 🔒
    
    # 输入输出 | Input/Output 📂
    input_paths: List[str] = field(default_factory=list)
    kmer_source_paths: List[str] = field(default_factory=list)
    query_target_paths: List[str] = field(default_factory=list)
    output_dir: str = "kmer_analysis_results"
    temp_dir: str = "/tmp/kmer_universal"
    
    # 角色分配 | Role Assignment 🏷️
    assignment_strategy: AssignmentStrategy = AssignmentStrategy.INTELLIGENT
    size_threshold_gb: float = 1.0   # 大于此值优先作为查询目标 ⚖️
    
    # 分析参数 | Analysis Parameters 🔬
    window_sizes: List[int] = field(default_factory=lambda: [500000])
    include_positions: bool = True   # 是否包含位置信息 📍
    output_formats: List[str] = field(default_factory=lambda: ["fasta", "csv", "txt"])
    
    # 性能优化 | Performance Optimization 🚀
    chunk_size_gb: float = 2.0      # 文件分片大小 🍰
    enable_compression: bool = True  # 启用压缩 📦
    cache_size_gb: float = 50.0     # 缓存大小 💨
    
    # 高级选项 | Advanced Options 🕹️
    keep_intermediate: bool = False  # 保留中间文件 🖇️
    verbose: bool = False           # 详细输出 🗣️
    dry_run: bool = False           # 干运行模式 💨
    resume: bool = False            # 断点续传 ⏯️
    
    def __post_init__(self):
        """配置后处理 ✅"""
        # 创建输出目录
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        Path(self.temp_dir).mkdir(parents=True, exist_ok=True)
        
        # 验证参数
        self.validate()
    
    def validate(self):
        """验证配置参数 ✔️"""
        if self.kmer_size < 1 or self.kmer_size > 256:
            raise ValueError(f"k-mer size must be between 1-256, got {self.kmer_size}")
        
        if self.threads < 1:
            raise ValueError(f"threads must be >= 1, got {self.threads}")
        
        if self.memory_gb < 1:
            raise ValueError(f"memory must be >= 1GB, got {self.memory_gb}")
        
        if self.signature_length not in range(5, 12):
            raise ValueError(f"signature length must be 5-11, got {self.signature_length}")
    
    @classmethod
    def from_yaml(cls, config_file: str) -> 'KmerConfig':
        """从YAML文件加载配置 📥📄"""
        with open(config_file, 'r', encoding='utf-8') as f:
            config_data = yaml.safe_load(f)
        
        return cls(**config_data)
    
    def to_yaml(self, config_file: str):
        """保存配置到YAML文件 📤📄"""
        config_dict = {
            'kmer_size': self.kmer_size,
            'threads': self.threads,
            'memory_gb': self.memory_gb,
            'min_count': self.min_count,
            'max_count': self.max_count,
            'canonical_form': self.canonical_form,
            'signature_length': self.signature_length,
            'ram_only_mode': self.ram_only_mode,
            'strict_memory': self.strict_memory,
            'input_paths': self.input_paths,
            'kmer_source_paths': self.kmer_source_paths,
            'query_target_paths': self.query_target_paths,
            'output_dir': self.output_dir,
            'temp_dir': self.temp_dir,
            'assignment_strategy': self.assignment_strategy.value,
            'size_threshold_gb': self.size_threshold_gb,
            'window_sizes': self.window_sizes,
            'include_positions': self.include_positions,
            'output_formats': self.output_formats,
            'chunk_size_gb': self.chunk_size_gb,
            'enable_compression': self.enable_compression,
            'cache_size_gb': self.cache_size_gb,
            'keep_intermediate': self.keep_intermediate,
            'verbose': self.verbose,
            'dry_run': self.dry_run,
            'resume': self.resume
        }
        
        with open(config_file, 'w', encoding='utf-8') as f:
            yaml.dump(config_dict, f, default_flow_style=False, allow_unicode=True)
    
    def get_kmc_params(self) -> Dict[str, Any]:
        """获取KMC参数字典 ⚙️"""
        return {
            'k': self.kmer_size,
            'm': self.memory_gb,
            'ci': self.min_count,
            'cx': self.max_count,
            'p': self.signature_length,
            't': self.threads,
            'b': not self.canonical_form,
            'r': self.ram_only_mode,
            'sm': self.strict_memory,
            'v': self.verbose
        }