"""
NeedleIdentity配置类|Needle Identity Configuration Class
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .utils import expand_path, get_tool_path


@dataclass
class NeedleIdentityConfig:
    """needle两两identity计算配置类|Needle Pairwise Identity Config Class"""

    # 必需参数|Required parameters
    input_file: str
    output_dir: str = "./output"

    # 可选参数|Optional parameters
    needle_path: Optional[str] = None
    threads: int = 12
    gapopen: float = 10.0
    gapextend: float = 0.5

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有~路径(关键,§11)|Expand all ~ paths (critical, §11)
        self.input_file = expand_path(self.input_file)
        self.output_dir = expand_path(self.output_dir)

        if self.needle_path:
            self.needle_path = expand_path(self.needle_path)
        else:
            # 优先级: 环境变量NEEDLE_PATH > 配置文件 > 默认值|Priority: NEEDLE_PATH env > config > default
            self.needle_path = get_tool_path(
                'needle',
                '~/miniforge3/envs/protein/bin/needle',
                'NEEDLE_PATH',
            )

        # 输出前缀|Output prefix
        self.output_prefix = Path(self.input_file).stem

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not self.input_file:
            errors.append("输入文件不能为空|Input file cannot be empty")
        elif not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")

        if not self.needle_path or not os.path.exists(self.needle_path):
            errors.append(f"needle不存在|needle not found: {self.needle_path}")

        if not isinstance(self.threads, int) or self.threads < 1:
            errors.append(f"线程数必须为正整数|Threads must be a positive integer: {self.threads}")

        if self.gapopen < 0:
            errors.append(f"gapopen必须非负|gapopen must be non-negative: {self.gapopen}")

        if self.gapextend < 0:
            errors.append(f"gapextend必须非负|gapextend must be non-negative: {self.gapextend}")

        if errors:
            raise ValueError("\n".join(errors))

        os.makedirs(self.output_dir, exist_ok=True)
        return True

    def __repr__(self):
        return (
            f"NeedleIdentityConfig(\n"
            f"  input_file={self.input_file!r},\n"
            f"  output_dir={self.output_dir!r},\n"
            f"  needle_path={self.needle_path!r},\n"
            f"  threads={self.threads!r},\n"
            f"  gapopen={self.gapopen!r},\n"
            f"  gapextend={self.gapextend!r}\n"
            f")"
        )
