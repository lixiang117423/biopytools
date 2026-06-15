"""插入检测模块配置类|Insert detection module configuration"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional
import os


@dataclass
class InsertDetectionConfig:
    """插入检测配置类|Insert detection configuration class"""

    # 必需参数|Required parameters
    genome: str
    insert_sequence: str
    fastq_dir: str
    output_dir: str

    # 可选参数|Optional parameters
    threads: int = 12
    skip_existing: bool = True  # 跳过已存在的文件|Skip existing files

    # 检测参数|Detection parameters
    min_clip: int = 20              # 最小soft-clip长度
    min_support: int = 5           # 最小支持reads数
    border_window: int = 500       # 边界窗口大小（bp）
    min_mapq: int = 10             # 最小比对质量
    score_threshold: int = 1000    # 得分阈值

    # 文件模式|File patterns（默认匹配fastp输出|Default matches fastp output）
    read1_suffix: str = "_1.clean.fq.gz"  # R1后缀（包含扩展名）
    read2_suffix: str = "_2.clean.fq.gz"  # R2后缀（包含扩展名）

    # 工具路径|Tool paths
    bowtie2_path: str = "bowtie2"
    samtools_path: str = "samtools"
    minimap2_path: str = "minimap2"  # 可选，用于长序列

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件是否存在|Check if input files exist
        if not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome}")

        if not os.path.exists(self.insert_sequence):
            errors.append(f"插入序列文件不存在|Insert sequence file not found: {self.insert_sequence}")

        if not os.path.exists(self.fastq_dir):
            errors.append(f"FASTQ目录不存在|FASTQ directory not found: {self.fastq_dir}")

        # 检查参数合法性|Check parameter validity
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive")

        if self.min_clip < 0:
            errors.append(f"最小clip长度必须>=0|min_clip must be >= 0")

        if self.min_support < 1:
            errors.append(f"最小支持reads数必须>=1|min_support must be >= 1")

        if errors:
            raise ValueError("\n".join(errors))

        return True
