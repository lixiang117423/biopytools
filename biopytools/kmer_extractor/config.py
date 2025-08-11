"""
🔧 K-mer提取配置管理模块 | K-mer Extraction Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union

@dataclass
class KmerConfig:
    """📋 K-mer提取配置类 | K-mer Extraction Configuration Class"""
    
    # 输入文件 | Input files
    input_files: Union[str, List[str]]
    output_dir: str = './kmer_output'
    
    # K-mer参数 | K-mer parameters
    kmer_length: int = 51  # 🧬 K-mer长度 | K-mer length
    
    # 性能参数 | Performance parameters
    threads: int = 88  # 🚀 线程数 | Number of threads
    memory_gb: int = 880  # 💾 内存限制(GB) | Memory limit (GB)
    
    # 文件类型和匹配 | File type and matching
    file_type: Optional[str] = None  # 'fasta' or 'fastq', auto-detect if None
    fastq_pattern: Optional[str] = None  # 📝 FASTQ文件匹配模式 | FASTQ file matching pattern
    
    # 处理选项 | Processing options
    canonical: bool = True  # 🔄 使用canonical k-mer
    compress_output: bool = True  # 🗜️ 压缩输出文件
    output_bed: bool = False  # 📋 是否输出BED文件 | Whether to output BED file
    keep_binary: bool = True  # 🗃️ 是否保留二进制文件 | Whether to keep binary files
    
    # 工具路径 | Tool paths
    unikmer_path: str = 'unikmer'
    jellyfish_path: str = 'jellyfish'
    
    # Jellyfish特有参数 | Jellyfish-specific parameters
    jellyfish_hash_size: str = '10000M'  # 🗂️ Jellyfish哈希表大小 | Jellyfish hash table size
    
    # 内部属性 | Internal attributes
    base_name: str = 'kmers'
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        if isinstance(self.input_files, str):
            self.input_files = [self.input_files]
        
        self.input_files = [os.path.normpath(os.path.abspath(f)) for f in self.input_files]
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 验证k-mer长度 | Validate k-mer length
        if self.kmer_length > 64:
            raise ValueError(f"❌ unikmer最大支持64-mer | unikmer supports maximum 64-mer: {self.kmer_length}")
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查输入文件/目录 | Check input files/directories
        for input_path in self.input_files:
            if not os.path.exists(input_path):
                errors.append(f"❌ 输入路径不存在 | Input path does not exist: {input_path}")
        
        # 检查参数范围 | Check parameter ranges
        if self.kmer_length <= 0 or self.kmer_length > 64:
            errors.append(f"❌ K-mer长度必须在1-64之间 | K-mer length must be between 1-64: {self.kmer_length}")
        
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Number of threads must be positive: {self.threads}")
        
        if self.memory_gb <= 0:
            errors.append(f"❌ 内存限制必须为正数 | Memory limit must be positive: {self.memory_gb}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
