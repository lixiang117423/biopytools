"""phylo_trim 配置管理模块|phylo_trim Configuration Management Module

聚合 mafft_fasttree + trimal 两模块的参数|Aggregates params for mafft_fasttree + trimal
"""

import os
from dataclasses import dataclass, field
from typing import Optional

from ..common.paths import expand_path, get_tool_path
from ..trimal.config import TrimalConfig


@dataclass
class PhyloTrimConfig:
    """phylo_trim 配置类|phylo_trim Configuration Class"""

    # 必需输入|Required input
    input_file: str

    # 输出|Output
    output_dir: str = "./phylo_trim_output"

    # trimal 开关|trimal toggle
    skip_trimal: bool = False

    # trimal 参数|trimal params
    trimal_method: str = "automated1"
    gt_threshold: float = 0.9
    cons_threshold: int = 80
    trimal_format: str = "keep"

    # mafft/fasttree 参数(透传 PhyloTreeBuilder)|mafft/fasttree params (passed to PhyloTreeBuilder)
    seq_type: Optional[str] = None        # 'protein'|'nucleotide',None=自动检测|auto-detect
    threads: int = 88
    mafft_params: str = "--auto"
    fasttree_params: str = ""
    mafft_path: str = "mafft"
    fasttree_path: str = "fasttree"

    # 其他|Other
    sample_name: Optional[str] = None
    log_file: Optional[str] = None
    verbose: bool = False

    # trimal 工具路径(支持 TRIMAL_PATH)|trimal tool path (overridable via TRIMAL_PATH)
    trimal_path: str = field(
        default_factory=lambda: get_tool_path(
            'trimal', '~/miniforge3/envs/trimal_v.1.5.0/bin/trimal', 'TRIMAL_PATH'
        )
    )

    # 复用 trimal 的合法枚举|reuse trimal's valid enums
    VALID_TRIMAL_METHODS = TrimalConfig.VALID_METHODS
    VALID_TRIMAL_FORMATS = TrimalConfig.VALID_FORMATS
    VALID_SEQ_TYPES = ('protein', 'nucleotide')

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.input_file = os.path.normpath(os.path.abspath(expand_path(self.input_file)))
        self.output_dir = os.path.normpath(os.path.abspath(expand_path(self.output_dir)))
        self.trimal_path = expand_path(self.trimal_path)
        self.mafft_path = expand_path(self.mafft_path)
        self.fasttree_path = expand_path(self.fasttree_path)

        if self.log_file is not None:
            self.log_file = os.path.normpath(os.path.abspath(expand_path(self.log_file)))

        if self.sample_name is None:
            self.sample_name = os.path.splitext(os.path.basename(self.input_file))[0]

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")

        if self.trimal_method not in self.VALID_TRIMAL_METHODS:
            errors.append(
                f"非法 trimal 方法|Invalid trimal method '{self.trimal_method}', "
                f"可选|valid: {', '.join(self.VALID_TRIMAL_METHODS)}"
            )

        if self.trimal_format not in self.VALID_TRIMAL_FORMATS:
            errors.append(
                f"非法 trimal 格式|Invalid trimal format '{self.trimal_format}', "
                f"可选|valid: {', '.join(self.VALID_TRIMAL_FORMATS)}"
            )

        if not (0.0 <= self.gt_threshold <= 1.0):
            errors.append(f"gt_threshold 越界|out of range [0,1]: {self.gt_threshold}")

        if not (0 <= self.cons_threshold <= 100):
            errors.append(f"cons_threshold 越界|out of range [0,100]: {self.cons_threshold}")

        if self.seq_type is not None and self.seq_type not in self.VALID_SEQ_TYPES:
            errors.append(
                f"非法 seq_type|Invalid seq_type '{self.seq_type}', "
                f"可选|valid: {', '.join(self.VALID_SEQ_TYPES)}"
            )

        if self.threads <= 0:
            errors.append(f"线程数必须为正|Threads must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
