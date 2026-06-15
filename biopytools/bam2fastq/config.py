"""
BAM to FASTQ转换配置管理模块|BAM to FASTQ Conversion Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class BAM2FASTQConfig:
    """BAM to FASTQ转换配置类|BAM to FASTQ Conversion Configuration Class"""

    # 必需参数|Required parameters
    input_dir: str  # 可以是文件或目录|Can be file or directory
    output_dir: str

    # 可选参数|Optional parameters
    threads: int = 64  # 每个BAM文件转换使用的线程数|Threads per BAM file conversion
    jobs: int = 1  # 并行处理的BAM文件数量|Number of parallel BAM file processing

    # bam2fastq配置|bam2fastq configuration
    bam2fastq_path: str = 'bam2fastq'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.input_path = Path(self.input_dir)
        self.output_path = Path(self.output_dir)

        # 创建输出目录|Create output directory
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.input_dir = str(self.input_path.absolute())
        self.output_dir = str(self.output_path.absolute())

        # 判断输入是文件还是目录（支持软链接）|Check if input is file or directory (supports symlinks)
        # 如果是软链接或后缀是.bam，判断为单文件|If symlink or .bam suffix, treat as single file
        if self.input_path.is_symlink() or self.input_path.suffix.lower() == '.bam':
            self.is_single_file = True
        else:
            self.is_single_file = self.input_path.is_file()

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入路径（可以是文件或目录，支持软链接）|Check input path (can be file or directory, supports symlinks)
        # 使用 lexists 检查软链接本身是否存在|Use lexists to check if symlink itself exists
        import os
        if not os.path.lexists(self.input_dir):
            errors.append(f"输入路径不存在|Input path does not exist: {self.input_dir}")

        # 如果是文件，检查是否为.bam文件|If file, check if it's a .bam file
        if self.is_single_file and not self.input_path.suffix.lower() == '.bam':
            errors.append(f"输入文件必须是.bam格式|Input file must be .bam format: {self.input_dir}")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        # 检查并行任务数|Check parallel job count
        if self.jobs <= 0:
            errors.append(f"并行任务数必须为正数|Job count must be positive: {self.jobs}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
