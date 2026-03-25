"""
AGP转表格配置管理模块|AGP to Table Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class AGPConfig:
    """AGP转表格配置类|AGP to Table Configuration Class"""

    # 必需文件|Required files
    agp_file: str  # AGP文件路径|AGP file path
    output_file: str  # 输出表格文件路径|Output table file path

    # 输出格式配置|Output format configuration
    format: str = 'txt'  # 输出格式 (txt/tsv/csv)|Output format
    add_statistics: bool = False  # 是否添加统计信息|Whether to add statistics
    group_by_scaffold: bool = True  # 是否按scaffold分组|Whether to group by scaffold
    add_headers: bool = True  # 是否添加表头|Whether to add headers

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.agp_file = os.path.normpath(os.path.abspath(self.agp_file))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))

        # 确保输出目录存在|Ensure output directory exists
        output_dir = os.path.dirname(self.output_file)
        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.agp_file):
            errors.append(f"AGP文件不存在|AGP file does not exist: {self.agp_file}")

        # 检查输出格式|Check output format
        valid_formats = ['txt', 'tsv', 'csv', 'xlsx']
        if self.format not in valid_formats:
            errors.append(f"无效的输出格式|Invalid output format: {self.format} (必须是|must be one of {valid_formats})")

        if errors:
            raise ValueError("\n".join(errors))

        return True
