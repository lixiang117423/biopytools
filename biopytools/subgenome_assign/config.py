"""
亚基因组归属配置模块|Subgenome Assignment Configuration Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from ..common.paths import expand_path, get_tool_path


@dataclass
class SubgenomeAssignConfig:
    """亚基因组归属配置类|Subgenome Assignment Configuration Class"""

    # 必需参数|Required parameters
    target: str                              # 目标多倍体基因组|Target polyploid genome
    parents: Dict[str, List[str]]            # {亲本名: [hap1.fa, hap2.fa, ...]}
    output_dir: str = './subgenome_assign_output'

    # 比对参数|Alignment parameters
    preset: str = 'asm10'        # minimap2 -x 预设: asm5/asm10/asm20|preset
    threads: int = 12
    minimap2_secondary: bool = False  # 是否输出次要比对|output secondary alignments

    # 归属判定参数|Assignment threshold
    min_conf: float = 0.65       # 置信度阈值，低于此值标记 LOW_CONFIDENCE|confidence threshold

    # FASTA 拆分|FASTA splitting
    split_fasta: bool = True     # 是否输出拆分后的 FASTA|output split FASTAs
    keep_unassigned: bool = True  # 是否输出未归属染色体的 FASTA|output unassigned FASTA

    # 工具路径|Tool paths
    minimap2_path: str = field(
        default_factory=lambda: get_tool_path(
            'minimap2', '~/miniforge3/envs/cphasing/bin/minimap2', 'MINIMAP2_PATH'
        )
    )
    samtools_path: str = field(
        default_factory=lambda: get_tool_path(
            'samtools', '~/.local/bin/samtools', 'SAMTOOLS_PATH'
        )
    )

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 关键：展开所有含~的路径|CRITICAL: Expand all ~ paths
        self.target = expand_path(self.target)
        self.output_dir = expand_path(self.output_dir)
        self.minimap2_path = expand_path(self.minimap2_path)
        self.samtools_path = expand_path(self.samtools_path)

        # 展开每个亲本的 hap 路径|Expand each parent's hap paths
        self.parents = {
            name: [expand_path(p) for p in haps]
            for name, haps in self.parents.items()
        }

        # 输出根目录|Output root
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 子目录|Subdirectories
        self.pipeline_info_dir = self.output_path / '00_pipeline_info'
        self.alignment_dir = self.output_path / '01_alignment'
        self.assignment_dir = self.output_path / '02_assignment'
        self.split_dir = self.output_path / '03_split_fastas'
        self.log_dir = self.output_path / '99_logs'

        for d in [self.pipeline_info_dir, self.alignment_dir,
                  self.assignment_dir, self.split_dir, self.log_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # 目标基因组名（去后缀）|Target genome name (no extension)
        self.target_name = Path(self.target).name

    def validate(self) -> List[str]:
        """验证配置参数，返回错误列表（空表示通过）|Validate config, return list of errors"""
        errors = []

        # 目标文件存在性|Target file existence
        if not os.path.exists(self.target):
            errors.append(f"目标基因组不存在|Target genome not found: {self.target}")

        # 亲本至少 2 个|At least 2 parents
        if len(self.parents) < 2:
            errors.append(
                f"亲本数量必须 >= 2|Need at least 2 parents, got {len(self.parents)}"
            )

        # 每个亲本的 hap 文件存在性|Each parent's hap files
        for name, haps in self.parents.items():
            if not haps:
                errors.append(f"亲本|Parent '{name}' 没有提供 hap 文件|has no hap files")
            for h in haps:
                if not os.path.exists(h):
                    errors.append(f"亲本|Parent '{name}' hap 不存在|not found: {h}")

        # 工具可执行性|Tools
        if not os.path.exists(self.minimap2_path):
            errors.append(f"minimap2 不存在|minimap2 not found: {self.minimap2_path}")
        if not os.path.exists(self.samtools_path):
            errors.append(f"samtools 不存在|samtools not found: {self.samtools_path}")

        # 参数范围|Parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正|Threads must be positive: {self.threads}")
        if not 0 < self.min_conf < 1:
            errors.append(f"min_conf 必须在 (0, 1)|min_conf must be in (0, 1): {self.min_conf}")
        if self.preset not in ('asm5', 'asm10', 'asm20', 'asm25'):
            errors.append(
                f"无效 preset|Invalid preset: {self.preset}; "
                f"valid: asm5/asm10/asm20/asm25"
            )

        return errors

    def get_minimap2_options(self) -> list:
        """生成minimap2参数列表|Generate minimap2 argument list"""
        args = ['-x', self.preset, '-t', str(self.threads)]
        if not self.minimap2_secondary:
            args.append('--secondary=no')
        return args

    def to_param_dict(self) -> dict:
        """导出参数字典（用于软件版本记录）|Export param dict for version log"""
        return {
            'target': self.target,
            'parents': {k: v for k, v in self.parents.items()},
            'preset': self.preset,
            'threads': self.threads,
            'min_conf': self.min_conf,
            'split_fasta': self.split_fasta,
            'keep_unassigned': self.keep_unassigned,
        }
