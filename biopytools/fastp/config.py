"""
FASTP质控配置管理模块|FASTP Quality Control Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class FastpConfig:
    """FASTP质控配置类|FASTP Quality Control Configuration Class"""

    # 必需参数|Required parameters
    input_dir: str  # 可以是目录或文件路径|Can be directory or file path
    output_dir: str

    # 软件配置|Software configuration
    fastp_path: str = "fastp"

    # 处理参数|Processing parameters
    threads: int = 12
    quality_threshold: int = 30
    min_length: int = 50
    unqualified_percent: int = 40
    n_base_limit: int = 10

    # 文件模式|File patterns
    read1_suffix: str = "_1.fq.gz"
    read2_suffix: str = "_2.fq.gz"
    single_end: bool = False  # 是否为单末端模式|Whether to use single-end mode

    # 日志配置|Logging configuration
    log_level: str = "INFO"
    quiet: bool = False
    verbose: int = 0

    # 执行控制|Execution control
    force: bool = False
    dry_run: bool = False

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.input_path = Path(self.input_dir)
        self.output_path = Path(self.output_dir)
        self.report_path = self.output_path / "fastp_reports"

        # 判断是否为单文件模式|Check if single file mode
        self.is_single_file = self.input_path.is_file()

        # 标准化路径|Normalize paths
        self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入路径（目录或文件）|Check input path (directory or file)
        if not self.input_path.exists():
            errors.append(f"输入路径不存在|Input path does not exist: {self.input_dir}")
        elif not (self.input_path.is_dir() or self.input_path.is_file()):
            errors.append(f"输入路径必须是目录或文件|Input path must be directory or file: {self.input_dir}")

        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive integer: {self.threads}")

        if self.quality_threshold < 0 or self.quality_threshold > 50:
            errors.append(f"质量阈值应在0-50之间|Quality threshold should be between 0-50: {self.quality_threshold}")

        if self.min_length <= 0:
            errors.append(f"最小长度必须为正整数|Minimum length must be positive integer: {self.min_length}")

        if self.unqualified_percent < 0 or self.unqualified_percent > 100:
            errors.append(f"不合格百分比应在0-100之间|Unqualified percentage should be between 0-100: {self.unqualified_percent}")

        if self.n_base_limit < 0:
            errors.append(f"N碱基限制不能为负数|N base limit cannot be negative: {self.n_base_limit}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
