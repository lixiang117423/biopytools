"""
Pi4Gene配置管理模块|Pi4Gene Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path, get_tool_path


@dataclass
class Pi4GeneConfig:
    """Pi4Gene配置类|Pi4Gene Configuration Class"""

    # 必需参数|Required parameters
    input_file: str
    id_file: str
    output_dir: str

    # 可选参数|Optional parameters
    threads: int = 12

    # 工具路径|Tool paths
    mafft_path: str = field(
        default_factory=lambda: get_tool_path(
            'mafft',
            '~/miniforge3/envs/mafft_v.7.525/bin/mafft',
            'MAFFT_PATH'
        )
    )

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.input_file = expand_path(self.input_file)
        self.id_file = expand_path(self.id_file)
        self.output_dir = expand_path(self.output_dir)
        self.mafft_path = expand_path(self.mafft_path)

        self.output_path = Path(self.output_dir).resolve()
        self.output_path.mkdir(parents=True, exist_ok=True)

        (self.output_path / '00_pipeline_info').mkdir(parents=True, exist_ok=True)
        (self.output_path / '01_mafft').mkdir(parents=True, exist_ok=True)
        (self.output_path / '99_logs').mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not os.path.exists(self.input_file):
            errors.append(f"序列文件不存在|Input file not found: {self.input_file}")

        if not os.path.exists(self.id_file):
            errors.append(f"分组ID文件不存在|ID file not found: {self.id_file}")

        if not os.path.exists(self.mafft_path):
            errors.append(f"mafft路径不存在|mafft path not found: {self.mafft_path}")

        if self.threads < 1:
            errors.append(f"线程数必须>=1|Threads must be >= 1: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
