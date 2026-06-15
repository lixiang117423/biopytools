"""
基因组挂载率统计配置管理模块|Genome Mount Rate Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass
class GenomeMountRateConfig:
    """基因组挂载率统计配置类|Genome Mount Rate Configuration Class"""

    # 必需参数|Required parameters
    fasta_file: str
    number: int

    # 可选参数|Optional parameters
    sort_by_length: bool = False

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize path
        self.fasta_file = os.path.normpath(os.path.abspath(self.fasta_file))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查FASTA文件是否存在|Check if FASTA file exists
        if not os.path.exists(self.fasta_file):
            errors.append(f"FASTA文件不存在|FASTA file does not exist: {self.fasta_file}")

        # 检查number参数|Check number parameter
        if self.number <= 0:
            errors.append(f"序列数量必须为正数|Number of sequences must be positive: {self.number}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
