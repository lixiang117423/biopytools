"""
VCF2Gene配置管理模块|VCF2Gene Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class VCF2GeneConfig:
    """VCF2Gene配置类|VCF2Gene Configuration Class"""

    # 必需参数|Required parameters
    vcf_file: str
    gff_file: str
    output_file: str

    # 可选参数|Optional parameters
    threads: int = 12

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))
        self.gff_file = os.path.normpath(os.path.abspath(self.gff_file))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))

        # 确保输出目录存在|Ensure output directory exists
        output_dir = Path(self.output_file).parent
        output_dir.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查VCF文件|Check VCF file
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在|VCF file does not exist: {self.vcf_file}")

        # 检查GFF文件|Check GFF file
        if not os.path.exists(self.gff_file):
            errors.append(f"GFF文件不存在|GFF file does not exist: {self.gff_file}")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
