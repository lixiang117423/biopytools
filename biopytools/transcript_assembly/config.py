"""
转录本从头组装配置管理模块|Transcript De Novo Assembly Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List


@dataclass
class TranscriptAssemblyConfig:
    """转录本从头组装配置类|Transcript De Novo Assembly Configuration Class"""

    # 必需参数|Required parameters
    genome_file: str
    input_dir: str
    output_dir: str

    # 处理参数|Processing parameters
    threads: int = 12
    fastq_pattern: str = '*_1.clean.fq.gz'
    sample_timeout: int = 43200  # 单个样本处理超时时间（秒），默认12小时|Sample processing timeout in seconds (default: 12 hours)

    # 步骤控制|Step control
    step: Optional[int] = None  # 1-6, None表示运行全部步骤|None means run all steps

    # 日志选项|Logging options
    log_file: Optional[str] = None
    log_level: str = "INFO"

    # 高级选项|Advanced options
    verbose: bool = False
    quiet: bool = False
    dry_run: bool = False
    force: bool = False

    # 内部属性|Internal attributes
    samples: List[dict] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.genome_file):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome_file}")

        # 检查输入目录|Check input directory
        if not os.path.isdir(self.input_dir):
            errors.append(f"输入目录不存在|Input directory not found: {self.input_dir}")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive integer: {self.threads}")

        # 检查步骤参数|Check step parameter
        if self.step is not None and self.step not in [1, 2, 3, 4, 5, 6]:
            errors.append(f"无效的步骤编号|Invalid step number: {self.step} (应为1-6|should be 1-6)")

        if errors:
            raise ValueError("\n".join(errors))

        return True
