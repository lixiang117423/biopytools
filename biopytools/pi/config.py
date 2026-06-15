"""
Pi计算配置管理模块|Pi Calculation Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path, get_tool_path


@dataclass
class PiConfig:
    """Pi计算配置类|Pi Calculation Configuration Class"""

    # 必需参数|Required parameters
    vcf_file: str
    pop_file: str
    output_dir: str
    genome_fai: str

    # 窗口参数|Window parameters
    window_size: Optional[int] = None  # None表示全基因组模式|None means genome-wide mode
    window_step: Optional[int] = None  # 步长，None则等于window_size（无重叠）|Step, None=window_size (no overlap)

    # 默认滑窗参数（全基因组模式时同时运行的滑窗计算）|Default windowed params (run alongside genome-wide)
    default_window_size: int = 100000   # 100kb
    default_window_step: int = 100000   # 步长等于窗口大小，无重叠

    # 质控参数|QC parameters
    maf: float = 0.0
    max_missing: float = 1.0

    # 线程数|Number of threads
    threads: int = 12

    # 输出控制|Output control
    keep_intermediate: bool = False

    # 工具路径 - 使用get_tool_path，支持~展开，禁止硬编码绝对路径
    # Tool paths - use get_tool_path, support ~ expansion, no hardcoded absolute paths
    vcftools_path: str = field(
        default_factory=lambda: get_tool_path(
            'vcftools',
            '~/miniforge3/envs/Population_genetics/bin/vcftools',
            'VCFTOOLS_PATH'
        )
    )

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有路径|Expand all paths
        self.vcf_file = expand_path(self.vcf_file)
        self.pop_file = expand_path(self.pop_file)
        self.output_dir = expand_path(self.output_dir)
        self.genome_fai = expand_path(self.genome_fai)
        self.vcftools_path = expand_path(self.vcftools_path)

        # 创建输出目录|Create output directories
        self.output_path = Path(self.output_dir).resolve()
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 创建子目录|Create subdirectories
        (self.output_path / '00_pipeline_info').mkdir(parents=True, exist_ok=True)
        (self.output_path / '01_vcftools').mkdir(parents=True, exist_ok=True)
        (self.output_path / '99_logs').mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在|VCF file does not exist: {self.vcf_file}")

        if not os.path.exists(self.pop_file):
            errors.append(f"群体文件不存在|Population file does not exist: {self.pop_file}")

        if not os.path.exists(self.genome_fai):
            errors.append(f"基因组fai文件不存在|Genome fai file does not exist: {self.genome_fai}")

        # 检查工具路径|Check tool paths
        if not os.path.exists(self.vcftools_path):
            errors.append(
                f"vcftools路径不存在|vcftools path does not exist: {self.vcftools_path}"
            )

        # 检查质控参数范围|Check QC parameter ranges
        if not 0 <= self.maf <= 1:
            errors.append(f"MAF阈值必须在0-1之间|MAF threshold must be between 0-1: {self.maf}")

        if not 0 <= self.max_missing <= 1:
            errors.append(
                f"最大缺失率必须在0-1之间|Max missing rate must be between 0-1: {self.max_missing}"
            )

        # 检查窗口大小|Check window size
        if self.window_size is not None and self.window_size < 1:
            errors.append(
                f"窗口大小必须>=1|Window size must be >= 1: {self.window_size}"
            )

        # 检查线程数|Check threads
        if self.threads < 1:
            errors.append(f"线程数必须>=1|Threads must be >= 1: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
