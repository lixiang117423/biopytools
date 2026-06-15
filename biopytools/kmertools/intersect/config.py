"""
配置管理模块|Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class KmerIntersectConfig:
    """Kmer交集配置类|Kmer Intersection Configuration Class"""

    # 输入文件|Input files
    kmer_matrix: str  # kmer矩阵文件|Kmer matrix file
    kmer_fasta: str   # 目标kmer的fasta文件|Target kmer fasta file

    # 输出文件|Output file
    output_file: str  # 输出文件|Output file

    # 处理参数|Processing parameters
    threads: int = 12  # 线程数|Number of threads
    chunk_size: int = 100000  # 每块行数|Lines per chunk for progress display

    # 查询选项|Query options
    use_reverse_complement: bool = True  # 是否使用反向互补查询|Whether to use reverse complement query

    # 输出选项|Output options
    keep_not_found: bool = True  # 是否保留未找到的kmer|Whether to keep kmers not found
    output_format: str = 'tsv'  # 输出格式: tsv或csv|Output format: tsv or csv
    window_size: int = 100000  # 窗口大小：用于统计每个窗口内的kmer数|Window size for kmer statistics per window

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.kmer_matrix = os.path.normpath(os.path.abspath(self.kmer_matrix))
        self.kmer_fasta = os.path.normpath(os.path.abspath(self.kmer_fasta))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))

        # 创建输出目录|Create output directory
        output_dir = os.path.dirname(self.output_file)
        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)

        # 验证输出格式|Validate output format
        self.output_format = self.output_format.lower()
        if self.output_format not in ['tsv', 'csv']:
            raise ValueError(f"输出格式必须是'tsv'或'csv'|Output format must be 'tsv' or 'csv': {self.output_format}")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input files
        if not os.path.exists(self.kmer_matrix):
            errors.append(f"kmer矩阵文件不存在|Kmer matrix file does not exist: {self.kmer_matrix}")

        if not os.path.exists(self.kmer_fasta):
            errors.append(f"kmer fasta文件不存在|Kmer fasta file does not exist: {self.kmer_fasta}")

        # 检查输出目录是否可写|Check output directory is writable
        output_dir = os.path.dirname(self.output_file)
        if output_dir and not os.access(output_dir, os.W_OK):
            errors.append(f"输出目录不可写|Output directory is not writable: {output_dir}")

        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if self.chunk_size <= 0:
            errors.append(f"块大小必须为正数|Chunk size must be positive: {self.chunk_size}")

        if self.window_size <= 0:
            errors.append(f"窗口大小必须为正数|Window size must be positive: {self.window_size}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
