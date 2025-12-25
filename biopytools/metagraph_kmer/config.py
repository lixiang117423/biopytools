"""
⚙️ K-mer分析配置管理模块 | K-mer Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class KmerConfig:
    """K-mer分析配置类 | K-mer Analysis Configuration Class"""
    
    # 输入文件 | Input files
    reference: str  # 参考序列文件或文件夹 | Reference sequence file or folder
    query: str  # 查询文件 | Query file
    output_dir: str = './kmer_output'
    
    # K-mer参数 | K-mer parameters
    kmer_length: int = 31
    
    # 性能参数 | Performance parameters
    threads: int = 88
    memory_gb: int = 64  # KMC内存限制 | KMC memory limit
    
    # 高级选项 | Advanced options
    canonical: bool = True  # 使用canonical k-mer | Use canonical k-mer
    min_count: int = 1  # 最小k-mer计数 | Minimum k-mer count
    
    # 工具路径 | Tool paths
    metagraph_path: str = 'metagraph'
    kmc_path: str = 'kmc'
    kmc_tools_path: str = 'kmc_tools'
    
    # 内部属性 | Internal attributes
    base_name: str = 'reference'
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.reference = os.path.normpath(os.path.abspath(self.reference))
        self.query = os.path.normpath(os.path.abspath(self.query))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 创建临时目录 | Create temporary directory
        self.tmp_dir = self.output_path / "tmp"
        self.tmp_dir.mkdir(parents=True, exist_ok=True)
        
        # MetaGraph文件路径 | MetaGraph file paths
        self.dbg_file = self.output_path / f"{self.base_name}.dbg"
        self.anno_file = self.output_path / f"{self.base_name}.column.annodbg"
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查参考文件/文件夹 | Check reference file/folder
        if not os.path.exists(self.reference):
            errors.append(f"❌ 参考文件/文件夹不存在 | Reference file/folder does not exist: {self.reference}")
        
        # 检查查询文件 | Check query file
        if not os.path.exists(self.query):
            errors.append(f"❌ 查询文件不存在 | Query file does not exist: {self.query}")
        
        # 检查参数范围 | Check parameter ranges
        if self.kmer_length < 1 or self.kmer_length > 255:
            errors.append(f"❌ k-mer长度必须在1-255之间 | k-mer length must be between 1-255: {self.kmer_length}")
        
        if self.threads < 1:
            errors.append(f"❌ 线程数必须为正整数 | Thread count must be positive: {self.threads}")
        
        if self.memory_gb < 1:
            errors.append(f"❌ 内存限制必须为正整数 | Memory limit must be positive: {self.memory_gb}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
