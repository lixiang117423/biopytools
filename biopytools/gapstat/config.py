"""
基因组Gap统计配置管理模块|Genome Gap Statistics Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class GapStatConfig:
    """基因组Gap统计配置类|Genome Gap Statistics Configuration Class"""

    # 必需文件|Required files
    fasta_file: str  # 输入FASTA文件

    # 输出文件|Output file
    output_file: Optional[str] = None  # 输出文件路径（可选，默认输出到终端）

    # 过滤参数|Filter parameters
    min_n: int = 1  # 最小连续N数量

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.fasta_file = os.path.normpath(os.path.abspath(self.fasta_file))

        # 如果指定了输出文件，确保输出目录存在|Ensure output directory exists if specified
        if self.output_file:
            self.output_file = os.path.normpath(os.path.abspath(self.output_file))
            output_dir = os.path.dirname(self.output_file)
            if output_dir:
                Path(output_dir).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.fasta_file):
            errors.append(f"FASTA文件不存在|FASTA file does not exist: {self.fasta_file}")

        # 检查min_n参数|Check min_n parameter
        if self.min_n < 1:
            errors.append(f"min_n必须大于0|min_n must be greater than 0: {self.min_n}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
