"""
Kraken2配置管理模块|Kraken2 Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path


@dataclass
class Kraken2Config:
    """Kraken2分析配置类|Kraken2 Analysis Configuration Class"""

    # 必需参数|Required parameters
    input_dir: str
    db_path: str

    # 可选参数|Optional parameters
    output_dir: str = './kraken2_output'
    threads: int = 12
    read_len: int = 150
    confidence: float = 0.0
    bracken_level: str = 'S'
    bracken_threshold: int = 10
    run_bracken: bool = True
    r1_suffix: str = '_1.clean.fq.gz'
    r2_suffix: str = '_2.clean.fq.gz'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.input_dir = expand_path(self.input_dir)
        self.db_path = expand_path(self.db_path)
        self.output_dir = expand_path(self.output_dir)

        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not os.path.isdir(self.input_dir):
            errors.append(
                f"输入目录不存在|Input directory not found: {self.input_dir}"
            )

        if not os.path.isdir(self.db_path):
            errors.append(
                f"Kraken2数据库目录不存在|Kraken2 database directory not found: {self.db_path}"
            )

        if self.threads <= 0:
            errors.append(
                f"线程数必须为正数|Thread count must be positive: {self.threads}"
            )

        if self.bracken_level not in ('D', 'P', 'C', 'O', 'F', 'G', 'S', 'S1'):
            errors.append(
                f"无效的Bracken分类级别|Invalid Bracken level: {self.bracken_level} "
                f"(应为|should be D/P/C/O/F/G/S/S1)"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return True
