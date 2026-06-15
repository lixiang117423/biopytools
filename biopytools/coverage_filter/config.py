"""
覆盖度过滤配置管理模块|Coverage Filter Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass
class CoverageFilterConfig:
    """覆盖度过滤配置类|Coverage Filter Configuration Class"""

    # 必需文件|Required files
    bam_file: str
    fasta_file: str
    output_prefix: str
    output_dir: str = "."

    # 线程数|Number of threads
    threads: int = 12

    # 高质量覆盖度阈值（默认90%）|High quality coverage threshold (default: 90%)
    high_coverage: float = 90.0

    # 中等质量最小覆盖度（默认50%）|Medium quality minimum coverage (default: 50%)
    medium_cov_min: float = 50.0

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.bam_file = os.path.normpath(os.path.abspath(self.bam_file))
        self.fasta_file = os.path.normpath(os.path.abspath(self.fasta_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 创建输出目录|Create output directory
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        # 计算衍生阈值|Calculate derived thresholds
        self.high_coverage_min = self.high_coverage
        self.medium_coverage_min = self.medium_cov_min
        self.medium_coverage_max = self.high_coverage
        self.low_coverage_max = self.medium_cov_min

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.bam_file):
            errors.append(f"BAM文件不存在|BAM file not found: {self.bam_file}")

        if not os.path.exists(self.fasta_file):
            errors.append(f"FASTA文件不存在|FASTA file not found: {self.fasta_file}")

        # 检查阈值合理性|Check threshold validity
        if not (0 < self.high_coverage <= 100):
            errors.append(f"高质量覆盖度阈值必须在0-100之间|High coverage threshold must be 0-100")

        if not (0 <= self.medium_cov_min < self.high_coverage):
            errors.append(f"中等质量最小值必须小于高质量阈值|Medium coverage minimum must be less than high coverage threshold")

        if errors:
            raise ValueError("\n".join(errors))

        return True
