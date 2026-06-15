"""
BAM统计分析配置模块|BAM Statistics Analysis Configuration Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

from ..common.paths import get_tool_path, expand_path


@dataclass
class BAMStatsConfig:
    """BAM统计分析配置类|BAM Statistics Analysis Configuration Class"""

    # 输入|Input
    input_path: str = ''
    reference_file: Optional[str] = None
    bed_file: Optional[str] = None

    # 输出|Output
    output_file: str = 'bam_stats.summary.tsv'

    # 分析参数|Analysis parameters
    min_mapq: int = 20
    min_base_quality: int = 20
    max_insert_size: int = 1000
    coverage_bins: int = 100
    window_size: int = 1_000_000
    step_size: int = 100_000

    # 性能设置|Performance settings
    threads: int = 12
    max_workers: int = 16

    # 模块开关|Module switches
    skip_alignment: bool = False
    skip_coverage: bool = False
    skip_sequence: bool = False
    skip_insert: bool = False
    skip_duplicate: bool = False
    skip_variation: bool = False

    # 工具路径|Tool paths
    samtools_path: str = 'samtools'
    bedtools_path: str = 'bedtools'

    # 内部状态|Internal state
    bam_files: List[str] = field(default_factory=list)
    output_path: Path = field(default_factory=Path)
    output_dir: str = field(default_factory=str)
    prefix: str = field(default_factory=str)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.input_path = expand_path(self.input_path)
        if self.reference_file:
            self.reference_file = expand_path(self.reference_file)
        if self.bed_file:
            self.bed_file = expand_path(self.bed_file)

        self.samtools_path = get_tool_path(
            'samtools', self.samtools_path, 'SAMTOOLS_PATH'
        )
        self.bedtools_path = get_tool_path(
            'bedtools', self.bedtools_path, 'BEDTOOLS_PATH'
        )

        self.output_file = expand_path(self.output_file)
        self.output_path = Path(self.output_file)
        self.output_dir = str(self.output_path.parent)
        self.prefix = self.output_path.stem
        self.output_path.parent.mkdir(parents=True, exist_ok=True)

        self._discover_bam_files()

    def _discover_bam_files(self):
        """发现BAM文件|Discover BAM files"""
        if not self.input_path:
            return

        input_path = Path(self.input_path)
        if input_path.is_file():
            if input_path.suffix.lower() in ['.bam', '.sam']:
                self.bam_files = [str(input_path)]
        elif input_path.is_dir():
            bam_files = []
            for pattern in ['*.bam', '*.BAM']:
                bam_files.extend(input_path.glob(pattern))
            self.bam_files = sorted(str(f) for f in bam_files)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not self.input_path:
            errors.append("必须指定输入路径|Input path must be specified")
        elif not Path(self.input_path).exists():
            errors.append(
                f"输入路径不存在|Input path does not exist: {self.input_path}"
            )

        if not self.bam_files:
            errors.append(
                f"未找到BAM文件|No BAM files found in: {self.input_path}"
            )

        if self.reference_file and not Path(self.reference_file).exists():
            errors.append(
                f"参考基因组不存在|Reference file not found: {self.reference_file}"
            )

        if self.bed_file and not Path(self.bed_file).exists():
            errors.append(
                f"BED文件不存在|BED file not found: {self.bed_file}"
            )

        if not 0 <= self.min_mapq <= 60:
            errors.append(
                f"MAPQ阈值须在0-60之间|MAPQ must be 0-60: {self.min_mapq}"
            )

        if self.threads <= 0:
            errors.append(
                f"线程数必须为正整数|Threads must be positive: {self.threads}"
            )

        if self.max_workers <= 0:
            errors.append(
                f"并行样本数必须为正整数|max_workers must be positive: {self.max_workers}"
            )

        if self.window_size <= 0:
            errors.append(
                f"窗口大小必须为正整数|Window size must be positive"
            )

        if self.step_size <= 0:
            errors.append(
                f"步长必须为正整数|Step size must be positive"
            )

        if self.step_size > self.window_size:
            errors.append(
                f"步长不能大于窗口|Step size > window size: "
                f"{self.step_size} > {self.window_size}"
            )

        if errors:
            raise ValueError('\n'.join(errors))
