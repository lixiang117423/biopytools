"""
FOF生成配置管理模块|FOF Generation Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from ..common.paths import expand_path


@dataclass
class FofConfig:
    """FOF生成配置类|FOF Generation Configuration Class"""

    # 必需参数|Required parameters
    input_path: str
    output_file: str

    # 可选参数|Optional parameters
    suffixes: List[str] = field(default_factory=list)  # 文件后缀过滤列表|File suffix filter list
    recursive: bool = False  # 是否递归扫描子目录|Whether to scan subdirectories recursively

    # 内部状态|Internal state
    is_single_file: bool = False  # 输入是否为单个文件|Whether input is a single file

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有~路径|Expand all ~ paths
        self.input_path = expand_path(self.input_path)
        self.output_file = expand_path(self.output_file)

        # 转换为Path对象|Convert to Path objects
        self.input_path_obj = Path(self.input_path)
        self.output_file_obj = Path(self.output_file)

        # 判断输入是文件还是目录（支持软链接）|Check if input is file or directory
        if not os.path.lexists(self.input_path):
            self.is_single_file = False
        elif self.input_path_obj.is_symlink():
            self.is_single_file = self.input_path_obj.resolve().is_file()
        else:
            self.is_single_file = self.input_path_obj.is_file()

        # 确保输出目录存在|Ensure output directory exists
        self.output_file_obj.parent.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入路径|Check input path
        if not os.path.lexists(self.input_path):
            errors.append(f"输入路径不存在|Input path does not exist: {self.input_path}")

        # 检查后缀格式|Check suffix format
        for suffix in self.suffixes:
            if not suffix.startswith('.'):
                errors.append(f"后缀必须以.开头|Suffix must start with '.': {suffix}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
