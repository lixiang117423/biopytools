"""
⚙️ K-mer PAV分析配置管理模块 | K-mer PAV Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class PAVConfig:
    """🔧 PAV分析配置类 | PAV Analysis Configuration Class"""
    
    # 📁 输入文件 | Input files
    genome_file: str
    fastq_dir: str
    output_dir: str = './kmer_pav_output'
    
    # 🧮 K-mer参数 | K-mer parameters
    kmer_size: int = 51
    threads: int = 8
    
    # 📊 分析参数 | Analysis parameters
    canonical: bool = True  # 使用canonical k-mer | Use canonical k-mers
    sort_output: bool = True  # 排序输出 | Sort output
    
    # 📂 FASTQ文件模式 | FASTQ file pattern
    fastq_pattern: str = "*_1.fq.gz"  # FASTQ文件匹配模式 | FASTQ file matching pattern
    
    # 🪟 窗口分析参数 | Window analysis parameters
    window_size: int = 500000  # 窗口大小，默认500kb | Window size, default 500kb
    step_size: Optional[int] = None  # 步长，默认None表示不重叠 | Step size, None means non-overlapping
    
    # 🛠️ 工具路径 | Tool paths
    unikmer_path: str = 'unikmer'
    
    # 💾 内部属性 | Internal attributes
    base_name: str = 'kmer_pav'
    
    def __post_init__(self):
        """🔧 初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 📍 标准化路径 | Normalize paths
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.fastq_dir = os.path.normpath(os.path.abspath(self.fastq_dir))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 🪟 处理窗口分析参数 | Process window analysis parameters
        if self.step_size is None:
            self.step_size = self.window_size  # 默认不重叠 | Default non-overlapping
            self.overlapping = False
        else:
            self.overlapping = True
        
        # 📁 内部文件路径 | Internal file paths
        self.genome_kmers_file = self.output_path / f"genome_{self.kmer_size}mers.unik"
        self.genome_positions_file = self.output_path / f"genome_{self.kmer_size}mers_positions.bed"
        self.sample_list_file = self.output_path / "sample_list.txt"
        
        # 📂 创建子目录 | Create subdirectories
        self.fastq_kmers_dir = self.output_path / "fastq_kmers"
        self.presence_results_dir = self.output_path / "presence_results"
        self.window_results_dir = self.output_path / "window_analysis"
        self.logs_dir = self.output_path / "logs"
        
        for subdir in [self.fastq_kmers_dir, self.presence_results_dir, self.window_results_dir, self.logs_dir]:
            subdir.mkdir(parents=True, exist_ok=True)
    
    def validate(self):
        """🔍 验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 🧬 检查基因组文件 | Check genome file
        if not os.path.exists(self.genome_file):
            errors.append(f"❌ 基因组文件不存在 | Genome file does not exist: {self.genome_file}")
        
        # 📂 检查fastq目录 | Check fastq directory
        if not os.path.exists(self.fastq_dir):
            errors.append(f"❌ FASTQ目录不存在 | FASTQ directory does not exist: {self.fastq_dir}")
        
        # 🔢 检查参数范围 | Check parameter ranges
        if self.kmer_size <= 0:
            errors.append(f"❌ K-mer长度必须为正整数 | K-mer size must be positive: {self.kmer_size}")
        
        if self.threads <= 0:
            errors.append(f"❌ 线程数必须为正整数 | Thread count must be positive: {self.threads}")
        
        if self.window_size <= 0:
            errors.append(f"❌ 窗口大小必须为正整数 | Window size must be positive: {self.window_size}")
        
        if self.step_size is not None and self.step_size <= 0:
            errors.append(f"❌ 步长必须为正整数 | Step size must be positive: {self.step_size}")
        
        # 🔍 检查fastq模式 | Check fastq pattern
        if '*' not in self.fastq_pattern:
            errors.append(f"❌ FASTQ模式必须包含通配符'*' | FASTQ pattern must contain wildcard '*': {self.fastq_pattern}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
