"""
染色体重命名配置管理模块|Chromosome Rename Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class ChrRenameConfig:
    """染色体重命名配置类|Chromosome Rename Configuration Class"""

    # 必需文件|Required files
    ref_fasta: str  # 参考基因组FASTA文件|Reference genome FASTA file
    query_fasta: str  # 待重命名的基因组FASTA文件|Query genome FASTA file to rename

    # 路径配置|Path configuration
    minimap2_path: str = 'minimap2'  # minimap2软件路径|minimap2 software path
    output_dir: str = './chr_rename_output'  # 输出目录|Output directory

    # 比对参数|Alignment parameters
    preset: str = 'asm5'  # minimap2预设模式|minimap2 preset mode (asm5/asm10/asm20)
    threads: int = 12  # 线程数|Number of threads

    # 过滤参数|Filtering parameters
    min_identity: float = 0.9  # 最小序列一致性阈值|Minimum identity threshold (0-1)
    min_alignment_length: int = 100000  # 最小比对长度(bp)|Minimum alignment length (bp)
    # 注意：覆盖度使用自动迭代策略（90% -> 20%），无需手动指定|Note: coverage uses auto iterative strategy (90% -> 20%), no manual specification needed

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.ref_fasta = os.path.normpath(os.path.abspath(self.ref_fasta))
        self.query_fasta = os.path.normpath(os.path.abspath(self.query_fasta))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 验证参数范围|Validate parameter ranges
        if not 0 < self.min_identity <= 1:
            raise ValueError(f"最小一致性必须在0-1之间|Min identity must be between 0-1: {self.min_identity}")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        required_files = [
            ('参考基因组文件|Reference genome file', self.ref_fasta),
            ('待重命名基因组文件|Query genome file', self.query_fasta),
        ]

        for file_desc, file_path in required_files:
            if not os.path.exists(file_path):
                errors.append(f"{file_desc}不存在|does not exist: {file_path}")

        # 检查minimap2是否可用|Check if minimap2 is available
        if self.minimap2_path != 'minimap2' and not os.path.exists(self.minimap2_path):
            errors.append(f"minimap2路径不存在|minimap2 path does not exist: {self.minimap2_path}")

        # 检查预设模式|Check preset mode
        valid_presets = ['asm5', 'asm10', 'asm20']
        if self.preset not in valid_presets:
            errors.append(f"无效的预设模式|Invalid preset mode: {self.preset} (必须是|must be one of {valid_presets})")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
