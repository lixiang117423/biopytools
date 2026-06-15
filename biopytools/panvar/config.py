"""泛基因组变异分析配置管理模块|Pan-genome Variant Analysis Configuration Management Module"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path, get_tool_path


@dataclass
class PanvarConfig:
    """泛基因组变异分析配置类|Pan-genome Variant Analysis Configuration Class"""

    # ============ 输入文件 ============
    input_file: str = ""                # GBZ图文件或VCF文件 (自动识别)
    ref_path: str = ""                  # vg deconstruct参考路径前缀 (如T2T, GBZ输入时必需)

    # ============ 输出 ============
    output_dir: str = "./panvar_output"

    # ============ 软件路径 ============
    vg_env: str = "vg_v.1.7.0"
    r_path: str = field(
        default_factory=lambda: get_tool_path(
            'Rscript',
            'Rscript',
            'R_PATH',
        )
    )

    # ============ 参数 ============
    threads: int = 12
    ref_size_mb: float = 0.0            # 参考基因组大小Mb; 0=自动推断
    n_permutations: int = 100           # 增长曲线随机置换次数

    # ============ 内部路径 (自动生成) ============
    input_type: str = field(default="", init=False)
    deconstruct_dir: Optional[str] = field(default=None, init=False)
    summarize_dir: Optional[str] = field(default=None, init=False)
    logs_dir: Optional[str] = field(default=None, init=False)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.input_file = expand_path(self.input_file) if self.input_file else ""
        self.output_dir = expand_path(self.output_dir) if self.output_dir else ""
        self.r_path = expand_path(self.r_path)

        # 自动识别输入类型|Auto-detect input type
        if self.input_file.endswith('.gbz'):
            self.input_type = 'gbz'
        else:
            self.input_type = 'vcf'

        # 创建输出目录结构|Create output directory structure
        out = Path(self.output_dir)
        self.deconstruct_dir = str(out / "01_deconstruct")
        self.summarize_dir = str(out / "02_summarize")
        self.logs_dir = str(out / "99_logs")

        for d in [self.deconstruct_dir, self.summarize_dir, self.logs_dir]:
            Path(d).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []
        if not self.input_file:
            errors.append("需要输入文件 -i/--input|Input file required")
        elif not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")
        if self.input_type == 'gbz' and not self.ref_path:
            errors.append("GBZ输入需要参考路径前缀 -P/--ref-path|GBZ input requires -P/--ref-path")
        if errors:
            raise ValueError("\n".join(errors))
        return True
