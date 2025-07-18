"""
K-mer分析配置管理模块 (双阶段设计) | K-mer Analysis Configuration Management Module (Two-Stage Design)
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class KmerConfig:
    """K-mer分析配置类 (双阶段设计) | K-mer Analysis Configuration Class (Two-Stage Design)"""
    
    # 双输入参数 | Dual input parameters
    database_input: str  # 数据库输入文件/目录 | Database input file/directory
    query_input: str     # 查询输入文件/目录 | Query input file/directory
    
    # 输出配置 | Output configuration
    output_prefix: str = "kmer_analysis"
    output_dir: str = "./kmer_output"
    
    # K-mer参数 | K-mer parameters
    kmer_size: int = 31
    min_count: int = 1  # 最小计数阈值 | Minimum count threshold
    max_count: int = 1000000  # 最大计数阈值 | Maximum count threshold
    reverse_complement: bool = False  # 是否包含反向互补 | Include reverse complement
    
    # 文件处理参数 | File processing parameters
    database_pattern: Optional[str] = None  # 数据库文件模式 | Database file pattern
    query_pattern: Optional[str] = None     # 查询文件模式 | Query file pattern
    threads: int = 8  # 线程数 | Number of threads
    
    # KMC特定参数 | KMC specific parameters
    kmc_memory_gb: int = 16  # KMC内存限制(GB) | KMC memory limit (GB)
    kmc_tmp_dir: str = "kmc_tmp"  # KMC临时目录 | KMC temporary directory
    keep_intermediate: bool = False  # 保留中间文件 | Keep intermediate files
    
    # 工具路径 | Tool paths
    kmc_path: str = "kmc"  # KMC可执行文件路径 | KMC executable path
    kmc_tools_path: str = "kmc_tools"  # kmc_tools可执行文件路径 | kmc_tools executable path
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 创建临时目录 | Create temporary directory
        self.tmp_path = self.output_path / self.kmc_tmp_dir
        self.tmp_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.database_input = os.path.normpath(os.path.abspath(self.database_input))
        self.query_input = os.path.normpath(os.path.abspath(self.query_input))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查输入路径 | Check input paths
        if not os.path.exists(self.database_input):
            errors.append(f"数据库输入路径不存在 | Database input path does not exist: {self.database_input}")
        
        if not os.path.exists(self.query_input):
            errors.append(f"查询输入路径不存在 | Query input path does not exist: {self.query_input}")
        
        # 检查k-mer大小 | Check k-mer size
        if self.kmer_size < 1 or self.kmer_size > 255:
            errors.append(f"K-mer大小必须在1-255之间 | K-mer size must be between 1-255: {self.kmer_size}")
        
        # 检查线程数 | Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Thread count must be positive integer: {self.threads}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
