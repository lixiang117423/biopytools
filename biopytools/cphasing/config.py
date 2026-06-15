"""
CPhasing配置管理模块|CPhasing Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List
from ..common.paths import expand_path, get_tool_path


@dataclass
class CPhasingConfig:
    """
    CPhasing配置类|CPhasing Configuration Class

    支持两种使用模式|Supports two modes:
    1. pipeline模式：完整封装的pipeline子命令|Fully wrapped pipeline subcommand
    2. 通用模式：透传任意子命令和参数|Pass-through any subcommand and arguments
    """

    # 子命令|Subcommand (默认pipeline，支持所有CPhasing子命令)
    subcommand: str = "pipeline"

    # === pipeline 专用参数 (pipeline subcommand parameters) ===
    fasta: Optional[str] = None       # 基因组FASTA文件|Genome FASTA file
    hic1: Optional[str] = None       # Hi-C R1 reads|Hi-C read1
    hic2: Optional[str] = None       # Hi-C R2 reads|Hi-C read2

    groups: str = "0"                # 分组数|Number of groups
    threads: int = 12                # 线程数|Number of threads
    mode: str = "phasing"            # 模式|Mode
    preset: str = "precision"        # 预设|Preset
    output_dir: str = "./cphasing_output"  # 输出目录|Output directory

    steps: Optional[str] = None      # 运行指定步骤|Run specified steps
    skip_steps: Optional[str] = None  # 跳过步骤|Skip steps

    hic_aligner: str = "_chromap"    # Hi-C比对器|Hi-C aligner
    hic_mapper_k: Optional[int] = None
    hic_mapper_w: Optional[int] = None
    mapping_quality: int = 0
    hcr: bool = False
    pattern: Optional[str] = None
    low_memory: bool = False

    # === 通用参数 (generic parameters) ===
    extra_args: List[str] = field(default_factory=list)  # 透传给CPhasing的额外参数

    # === CPhasing软件目录 ===
    cphasing_dir: str = field(
        default_factory=lambda: get_tool_path(
            'cphasing',
            '~/software/CPhasing/CPhasing-main',
            'CPHASING_DIR'
        )
    )

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.cphasing_dir = expand_path(self.cphasing_dir)
        self.output_dir = expand_path(self.output_dir)

        for attr in ('fasta', 'hic1', 'hic2'):
            val = getattr(self, attr)
            if val is not None:
                val = expand_path(val)
                if not os.path.isabs(val):
                    val = str(Path(val).resolve())
                setattr(self, attr, val)

        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # pipeline子命令需要必需文件|pipeline subcommand requires input files
        if self.subcommand == "pipeline":
            if not self.fasta:
                errors.append("pipeline模式需要--fasta参数|pipeline mode requires --fasta")
            if not self.hic1:
                errors.append("pipeline模式需要--hic1参数|pipeline mode requires --hic1")
            if not self.hic2:
                errors.append("pipeline模式需要--hic2参数|pipeline mode requires --hic2")

            if self.fasta and not os.path.exists(self.fasta):
                errors.append(f"基因组FASTA文件不存在|Genome FASTA not found: {self.fasta}")
            if self.hic1 and not os.path.exists(self.hic1):
                errors.append(f"Hi-C R1文件不存在|Hi-C R1 not found: {self.hic1}")
            if self.hic2 and not os.path.exists(self.hic2):
                errors.append(f"Hi-C R2文件不存在|Hi-C R2 not found: {self.hic2}")

            valid_modes = ['phasing', 'haploid', 'basal', 'basal_withprune']
            if self.mode not in valid_modes:
                errors.append(f"无效的模式|Invalid mode: {self.mode}")

            valid_presets = ['precision', 'sensitive', 'very-sensitive', 'nofilter']
            if self.preset not in valid_presets:
                errors.append(f"无效的预设|Invalid preset: {self.preset}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
