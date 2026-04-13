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
    # None 表示自动检测，支持自动识别 .fq.gz 和 .fastq.gz
    # None means auto-detect, supports both .fq.gz and .fastq.gz
    read1_suffix: Optional[str] = None
    read2_suffix: Optional[str] = None
    single_end: bool = False  # 是否为单末端模式|Whether to use single-end mode

    # 日志配置|Logging configuration
    log_level: str = "INFO"
    quiet: bool = False
    verbose: int = 0

    # 执行控制|Execution control
    force: bool = False
    dry_run: bool = False

    # SeqKit配对修复配置|SeqKit pair configuration
    enable_pair: bool = True  # 是否启用seqkit pair配对修复步骤（默认启用）|Whether to enable seqkit pair step (enabled by default)
    seqkit_path: str = "seqkit"  # seqkit可执行文件路径|seqkit executable path

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.input_path = Path(self.input_dir)
        self.output_path = Path(self.output_dir)
        self.report_path = self.output_path / "fastp_reports"

        # 判断是否为单文件模式|Check if single file mode
        self.is_single_file = self.input_path.is_file()

        # 自动检测后缀（如果用户未指定）|Auto-detect suffixes if not specified by user
        if self.read1_suffix is None or self.read2_suffix is None:
            self._auto_detect_suffixes()

    def _auto_detect_suffixes(self):
        """
        自动检测输入文件的后缀格式|Auto-detect input file suffix format

        支持 .fq.gz 和 .fastq.gz 格式的自动识别
        Supports both .fq.gz and .fastq.gz formats

        优先顺序|Priority order:
        1. _1.fq.gz / _2.fq.gz
        2. _1.fastq.gz / _2.fastq.gz
        """
        # 定义可能的后缀组合|Define possible suffix combinations
        possible_suffixes = [
            ("_1.fq.gz", "_2.fq.gz"),
            ("_1.fastq.gz", "_2.fastq.gz"),
        ]

        # 检测目录模式|Detect in directory mode
        if not self.is_single_file and self.input_path.is_dir():
            detected = False
            for read1_suf, read2_suf in possible_suffixes:
                # 尝试查找匹配 read1 后缀的文件|Try to find files matching read1 suffix
                test_files = list(self.input_path.glob(f"*{read1_suf}"))
                if test_files:
                    self.read1_suffix = read1_suf
                    self.read2_suffix = read2_suf
                    detected = True
                    break

            if detected and self.read1_suffix:
                import logging
                logger = logging.getLogger(__name__)
                logger.info(f"自动检测到文件后缀|Auto-detected file suffix: {self.read1_suffix}")
                return

        # 检测单文件模式|Detect in single file mode
        if self.is_single_file:
            filename = self.input_path.name
            for read1_suf, read2_suf in possible_suffixes:
                if filename.endswith(read1_suf):
                    self.read1_suffix = read1_suf
                    self.read2_suffix = read2_suf

                    import logging
                    logger = logging.getLogger(__name__)
                    logger.info(f"自动检测到文件后缀|Auto-detected file suffix: {self.read1_suffix}")
                    return

        # 如果都检测不到，使用默认值|If none detected, use default
        if self.read1_suffix is None:
            self.read1_suffix = "_1.fq.gz"
        if self.read2_suffix is None:
            self.read2_suffix = "_2.fq.gz"

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
