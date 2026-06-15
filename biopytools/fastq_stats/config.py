"""
FASTQ文件统计配置管理模块|FASTQ File Statistics Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class FastqStatsConfig:
    """FASTQ文件统计配置类|FASTQ File Statistics Configuration Class"""

    # 必需参数|Required parameters
    input_path: str
    output_file: str

    # 可选参数|Optional parameters
    pattern: Optional[str] = None
    threads: int = 12

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.input_path = os.path.normpath(os.path.abspath(self.input_path))
        self.output_file = os.path.normpath(os.path.abspath(self.output_file))

        # 获取输出目录|Get output directory
        self.output_dir = os.path.dirname(self.output_file)
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 验证输出文件格式|Validate output file format
        if not (self.output_file.endswith('.csv') or self.output_file.endswith('.xlsx')):
            raise ValueError(
                f"输出文件格式不支持|Unsupported output file format: {self.output_file}. "
                f"仅支持.csv和.xlsx格式|Only .csv and .xlsx formats are supported"
            )

        # 确定输出格式|Determine output format
        self.output_format = 'excel' if self.output_file.endswith('.xlsx') else 'csv'

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入路径|Check input path
        if not os.path.exists(self.input_path):
            errors.append(
                f"输入路径不存在|Input path does not exist: {self.input_path}"
            )

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(
                f"线程数必须为正数|Thread count must be positive: {self.threads}"
            )

        # 检查pattern格式|Check pattern format
        if self.pattern and '*' not in self.pattern:
            errors.append(
                f"模式必须包含*通配符|Pattern must contain * wildcard: {self.pattern}"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return True
