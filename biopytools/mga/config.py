"""MGA共识组装配置模块|MGA consensus assembly configuration module"""

import os
from dataclasses import dataclass
from pathlib import Path

from ..common.paths import expand_path, get_tool_path


@dataclass
class MGAConfig:
    """MGA配置类|MGA configuration class"""

    reads: str
    output_dir: str
    threads: int = 50
    mga_path: str = None
    conda_env: str = "mga"
    dry_run: bool = False

    def __post_init__(self):
        # mga_path默认走路径优先级(env var MGA_PATH > 配置文件 > 默认)|Default via path priority
        if self.mga_path is None:
            self.mga_path = get_tool_path(
                "mga", "~/software/MGA/consensusLJA/bin/MGA", "MGA_PATH"
            )
        # 展开所有~路径(关键)|Expand all ~ paths (critical)
        self.reads = expand_path(self.reads)
        self.output_dir = expand_path(self.output_dir)
        self.mga_path = expand_path(self.mga_path)

        # 输出目录与规范子目录|Output dir and standard subdirs
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.info_path = self.output_path / "00_pipeline_info"
        self.info_path.mkdir(parents=True, exist_ok=True)
        self.logs_path = self.output_path / "99_logs"
        self.logs_path.mkdir(parents=True, exist_ok=True)

        # 最终产物(断点续传判断)|Final product (checkpoint check)
        self.assembly_fasta = self.output_path / "5_polishing" / "assembly.fasta"

    def validate(self):
        """验证配置|Validate configuration"""
        errors = []
        if not os.path.exists(self.reads):
            errors.append(f"输入reads不存在|Input reads not found: {self.reads}")
        if not os.path.exists(self.mga_path):
            errors.append(f"MGA二进制不存在|MGA binary not found: {self.mga_path}")
        if self.threads <= 0:
            errors.append(f"线程数必须为正|Threads must be positive: {self.threads}")
        if errors:
            raise ValueError("\n".join(errors))
        return True
