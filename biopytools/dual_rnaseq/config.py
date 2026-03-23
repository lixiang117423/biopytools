"""
互作转录组配置管理模块|Dual RNA-seq Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List


@dataclass
class DualRNASeqConfig:
    """互作转录组配置类|Dual RNA-seq Configuration Class"""

    # 物种1配置|Species 1 configuration (如宿主|e.g., host)
    species1_name: str
    species1_genome: str
    species1_gtf: str

    # 物种2配置|Species 2 configuration (如病原体|e.g., pathogen)
    species2_name: str
    species2_genome: str
    species2_gtf: str

    # 通用参数|Common parameters
    input_path: str
    output_dir: str

    # 处理参数|Processing parameters
    threads: int = 12
    fastq_pattern: Optional[str] = None  # FASTQ文件命名模式|FASTQ naming pattern
    extract_fastq: bool = True  # 是否从BAM提取FASTQ|Whether to extract FASTQ from BAM (default: True)

    # 分类参数|Classification parameters
    min_mapq: int = 20  # 最小mapping quality|Minimum mapping quality
    unique_only: bool = True  # 仅保留唯一比对|Only keep unique mappings

    # 内部属性|Internal attributes
    samples: List[dict] = None  # 样本信息列表|Sample information list

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.species1_genome = os.path.normpath(os.path.abspath(self.species1_genome))
        self.species1_gtf = os.path.normpath(os.path.abspath(self.species1_gtf))
        self.species2_genome = os.path.normpath(os.path.abspath(self.species2_genome))
        self.species2_gtf = os.path.normpath(os.path.abspath(self.species2_gtf))
        self.input_path = os.path.normpath(os.path.abspath(self.input_path))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查物种1文件|Check species 1 files
        species1_files = [
            ('物种1基因组文件|Species 1 genome file', self.species1_genome),
            ('物种1 GTF文件|Species 1 GTF file', self.species1_gtf),
        ]

        for file_desc, file_path in species1_files:
            if not os.path.exists(file_path):
                errors.append(f"{file_desc}不存在|does not exist: {file_path}")

        # 检查物种2文件|Check species 2 files
        species2_files = [
            ('物种2基因组文件|Species 2 genome file', self.species2_genome),
            ('物种2 GTF文件|Species 2 GTF file', self.species2_gtf),
        ]

        for file_desc, file_path in species2_files:
            if not os.path.exists(file_path):
                errors.append(f"{file_desc}不存在|does not exist: {file_path}")

        # 检查输入路径|Check input path
        if not os.path.exists(self.input_path):
            errors.append(f"输入路径不存在|Input path does not exist: {self.input_path}")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive integer: {self.threads}")

        # 检查min_mapq参数|Check min_mapq parameter
        if self.min_mapq < 0 or self.min_mapq > 60:
            errors.append(f"min_mapq必须在0-60之间|min_mapq must be between 0-60: {self.min_mapq}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
