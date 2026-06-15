"""
Kmer矩阵比较配置管理模块|Kmer Matrix Comparison Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class KmerCompareConfig:
    """Kmer矩阵比较配置类|Kmer Matrix Comparison Configuration Class"""

    # 必需参数|Required parameters
    file1: str  # 第一个kmer矩阵文件|First kmer matrix file
    file2: str  # 第二个kmer矩阵文件|Second kmer matrix file
    output_prefix: str  # 输出文件前缀|Output file prefix

    # 可选参数|Optional parameters
    window_size: int = 100000  # 窗口大小（行数）|Window size in lines

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.file1 = os.path.normpath(os.path.abspath(self.file1))
        self.file2 = os.path.normpath(os.path.abspath(self.file2))
        self.output_prefix = os.path.normpath(os.path.abspath(self.output_prefix))

        # 创建输出目录|Create output directory
        output_dir = os.path.dirname(self.output_prefix)
        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)

        # 输出文件名|Output file names
        self.output_file1 = f"{self.output_prefix}_file1_stats.txt"
        self.output_file2 = f"{self.output_prefix}_file2_stats.txt"
        self.log_file = f"{self.output_prefix}.log"

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input files
        if not os.path.exists(self.file1):
            errors.append(f"文件1不存在|File1 does not exist: {self.file1}")

        if not os.path.exists(self.file2):
            errors.append(f"文件2不存在|File2 does not exist: {self.file2}")

        # 检查参数范围|Check parameter ranges
        if self.window_size <= 0:
            errors.append(f"窗口大小必须为正数|Window size must be positive: {self.window_size}")

        # 检查输出目录是否可写|Check output directory is writable
        output_dir = os.path.dirname(self.output_prefix)
        if output_dir and not os.access(output_dir, os.W_OK):
            errors.append(f"输出目录不可写|Output directory is not writable: {output_dir}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
