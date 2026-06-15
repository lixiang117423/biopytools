"""
FASTQ GC过滤配置管理模块|FASTQ GC Filter Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class FastqGcFilterConfig:
    """FASTQ GC过滤配置类|FASTQ GC Filter Configuration Class"""

    # 必需参数|Required parameters
    input_file: str
    output_file: str

    # GC含量过滤参数|GC content filtering parameters
    min_gc: float = 25.0
    max_gc: float = 100.0

    # 序列长度过滤参数|Sequence length filtering parameters
    min_length: int = 50  # 最短序列长度|Minimum sequence length
    max_length: Optional[int] = None  # 最长序列长度, None表示不限制|Maximum sequence length, None means no limit

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))

        # 确保输出目录存在|Ensure output directory exists
        output_dir = os.path.dirname(self.output_file)
        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件是否存在|Check if input file exists
        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file does not exist: {self.input_file}")

        # 检查GC含量范围|Check GC content range
        if not (0 <= self.min_gc <= 100):
            errors.append(f"最小GC含量必须在0-100之间|Min GC content must be between 0-100: {self.min_gc}")
        if not (0 <= self.max_gc <= 100):
            errors.append(f"最大GC含量必须在0-100之间|Max GC content must be between 0-100: {self.max_gc}")
        if self.min_gc > self.max_gc:
            errors.append(f"最小GC含量不能大于最大GC含量|Min GC content cannot be greater than max GC content: {self.min_gc} > {self.max_gc}")

        # 检查序列长度范围|Check sequence length range
        if self.min_length < 0:
            errors.append(f"最小序列长度不能为负数|Min sequence length cannot be negative: {self.min_length}")
        if self.max_length is not None and self.max_length < 0:
            errors.append(f"最大序列长度不能为负数|Max sequence length cannot be negative: {self.max_length}")
        if self.max_length is not None and self.min_length > self.max_length:
            errors.append(f"最小序列长度不能大于最大序列长度|Min length cannot be greater than max length: {self.min_length} > {self.max_length}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
