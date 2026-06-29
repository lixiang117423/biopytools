"""a-liner pipeline 配置模块|aliner pipeline configuration module"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

from ..common.paths import expand_path, get_tool_path
from .utils import parse_seq_spec


@dataclass
class AlinerConfig:
    """a-liner pipeline 配置类|aliner pipeline configuration class"""

    # 必需参数|required parameters
    ref_fasta: str
    query_fasta: str
    ref_seqs: List[str]                 # 如 ['chrZ', 'chrZ:1-30000000']
    query_seqs: List[str]

    # 输出参数|output parameters
    output_dir: str = "./aliner_output"
    out_prefix: str = "synteny"

    # 比对参数|alignment parameters
    preset: str = "asm5"                # minimap2 预设|minimap2 preset
    min_identity: int = 70              # a-liner identity 阈值|a-liner identity threshold
    min_alignment_len: int = 1000       # a-liner 最小比对长度|a-liner min alignment length
    threads: int = 12

    # 可视化参数|visualization parameters
    colormap: int = 5
    figure_size: Tuple[float, float] = (6, 0)
    extra_args: str = ""                # 透传给a-liner的兜底参数|pass-through args to a-liner

    # 工具路径|tool paths
    minimap2_path: str = field(
        default_factory=lambda: get_tool_path('minimap2', '~/miniforge3/envs/telocomp/bin/minimap2', 'MINIMAP2_PATH'))
    samtools_path: str = field(
        default_factory=lambda: get_tool_path('samtools', '~/miniforge3/envs/telocomp/bin/samtools', 'SAMTOOLS_PATH'))
    aliner_env: str = "a-liner"         # a-liner固定conda环境|a-liner fixed conda env

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 关键：展开所有~路径|CRITICAL: expand all ~ paths
        self.ref_fasta = expand_path(self.ref_fasta)
        self.query_fasta = expand_path(self.query_fasta)
        self.output_dir = expand_path(self.output_dir)
        self.minimap2_path = expand_path(self.minimap2_path)
        self.samtools_path = expand_path(self.samtools_path)

        # 解析序列规格（非法格式在此抛ValueError）|parse seq specs (invalid raises here)
        self.ref_specs = [parse_seq_spec(s) for s in self.ref_seqs]
        self.query_specs = [parse_seq_spec(s) for s in self.query_seqs]

        # 创建输出目录|create output dirs
        for sub in ['00_pipeline_info', '01_alignment', '02_aliner', '99_logs']:
            Path(self.output_dir, sub).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置|Validate configuration"""
        errors = []
        if not os.path.exists(self.ref_fasta):
            errors.append(f"参考基因组不存在|Reference genome not found: {self.ref_fasta}")
        if not os.path.exists(self.query_fasta):
            errors.append(f"查询基因组不存在|Query genome not found: {self.query_fasta}")
        if len(self.ref_seqs) != len(self.query_seqs):
            errors.append(
                f"ref-seqs与query-seqs数量不一致（配对模式要求等长）|"
                f"ref-seqs and query-seqs length mismatch (paired mode requires equal length): "
                f"{len(self.ref_seqs)} != {len(self.query_seqs)}")
        if self.threads <= 0:
            errors.append(f"线程数必须为正|Threads must be positive: {self.threads}")
        if self.preset not in ('asm5', 'asm10', 'asm20'):
            errors.append(f"不支持的preset|Unsupported preset: {self.preset}（asm5/asm10/asm20）")
        if errors:
            raise ValueError("\n".join(errors))
        return True
