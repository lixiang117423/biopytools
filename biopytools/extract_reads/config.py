"""
Reads提取配置管理模块|Reads Extraction Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class ExtractReadsConfig:
    """Reads提取配置类|Reads Extraction Configuration Class"""

    # 必需参数|Required parameters
    mapping_file: str
    fastq_file: str
    output_file: str

    # 可选参数|Optional parameters
    compress_output: bool = True

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.mapping_file = os.path.normpath(os.path.abspath(self.mapping_file))
        self.fastq_file = os.path.normpath(os.path.abspath(self.fastq_file))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))

        # 如果没有指定.gz后缀且需要压缩，自动添加|Add .gz extension if compress and not specified
        if self.compress_output and not self.output_file.endswith('.gz'):
            self.output_file = self.output_file + '.gz'

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.mapping_file):
            errors.append(f"对应关系表不存在|Mapping file not found: {self.mapping_file}")

        if not os.path.exists(self.fastq_file):
            errors.append(f"FASTQ文件不存在|FASTQ file not found: {self.fastq_file}")

        # 检查文件扩展名|Check file extensions
        if not self.fastq_file.endswith(('.fq', '.fq.gz', '.fastq', '.fastq.gz')):
            errors.append(f"FASTQ文件扩展名错误|Invalid FASTQ file extension: {self.fastq_file}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
