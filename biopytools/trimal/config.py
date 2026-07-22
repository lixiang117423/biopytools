"""trimal 配置管理模块|trimal Configuration Management Module"""

import os
from dataclasses import dataclass, field
from typing import Optional

from ..common.paths import expand_path, get_tool_path


@dataclass
class TrimalConfig:
    """trimal 配置类|trimal Configuration Class"""

    # 必需输入|Required input
    input_file: str

    # 输出|Output
    output_dir: str = "./trimal_output"

    # 修剪方法|Trimming method
    method: str = "automated1"
    gt_threshold: float = 0.9       # 仅 method=gt 生效,范围[0,1]|Only for method=gt, range [0,1]
    cons_threshold: int = 80        # 仅 method=cons 生效,范围[0,100]|Only for method=cons, range [0,100]

    # 输出格式|Output format
    output_format: str = "keep"     # keep=沿用输入格式|keep=input format

    # 附加输出|Extra outputs
    colnumbering: bool = False
    backtrans_file: Optional[str] = None
    complementary: bool = False

    # 其他|Other
    keep_header: bool = False
    sample_name: Optional[str] = None
    log_file: Optional[str] = None
    verbose: bool = False

    # 工具路径(支持 TRIMAL_PATH 环境变量覆盖)|Tool path (overridable via TRIMAL_PATH)
    trimal_path: str = field(
        default_factory=lambda: get_tool_path(
            'trimal',
            '~/miniforge3/envs/trimal_v.1.5.0/bin/trimal',
            'TRIMAL_PATH'
        )
    )

    VALID_METHODS = ('automated1', 'gappyout', 'strict', 'strictplus', 'gt', 'cons')
    VALID_FORMATS = ('keep', 'fasta', 'phylip', 'phylip_paml', 'clustal', 'nexus', 'nbrf', 'mega')

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开并标准化路径|Expand and normalize paths
        self.input_file = os.path.normpath(os.path.abspath(expand_path(self.input_file)))
        self.output_dir = os.path.normpath(os.path.abspath(expand_path(self.output_dir)))
        self.trimal_path = expand_path(self.trimal_path)

        if self.backtrans_file is not None:
            self.backtrans_file = os.path.normpath(os.path.abspath(expand_path(self.backtrans_file)))

        if self.log_file is not None:
            self.log_file = os.path.normpath(os.path.abspath(expand_path(self.log_file)))

        # 默认 sample_name 取输入文件名(去扩展名)|Default sample_name from input basename
        if self.sample_name is None:
            self.sample_name = os.path.splitext(os.path.basename(self.input_file))[0]

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 输入文件|Input file
        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")

        # 方法|Method
        if self.method not in self.VALID_METHODS:
            errors.append(
                f"非法方法|Invalid method '{self.method}', 可选|valid: {', '.join(self.VALID_METHODS)}"
            )

        # 格式|Format
        if self.output_format not in self.VALID_FORMATS:
            errors.append(
                f"非法格式|Invalid format '{self.output_format}', 可选|valid: {', '.join(self.VALID_FORMATS)}"
            )

        # 阈值范围|Threshold ranges
        if not (0.0 <= self.gt_threshold <= 1.0):
            errors.append(
                f"gt_threshold 越界|out of range [0,1]: {self.gt_threshold}"
            )

        if not (0 <= self.cons_threshold <= 100):
            errors.append(
                f"cons_threshold 越界|out of range [0,100]: {self.cons_threshold}"
            )

        # backtrans 文件|backtrans file
        if self.backtrans_file is not None and not os.path.exists(self.backtrans_file):
            errors.append(f"backtrans CDS 文件不存在|backtrans CDS file not found: {self.backtrans_file}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
