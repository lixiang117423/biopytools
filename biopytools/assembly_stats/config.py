"""
基因组装配统计配置管理模块|Genome Assembly Statistics Configuration Management Module
"""

import os
import glob
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


@dataclass
class AssemblyStatsConfig:
    """基因组装配统计配置类|Genome Assembly Statistics Configuration Class"""

    # 输入路径|Input path (file or directory)
    input_path: str

    # 输出配置|Output configuration
    output_dir: str = "./assembly_stats_output"
    output_prefix: str = "assembly_stats"

    # 过滤参数|Filter parameters
    min_length: int = 1
    file_pattern: str = "*.fa"  # 文件匹配模式|File pattern for directory

    # 输出格式|Output format
    grep_friendly: bool = False
    tab_delimited: bool = False
    no_header: bool = False
    vertical_format: bool = True  # 竖着格式|Vertical format (key-value pairs)

    # 报告格式|Report format
    generate_csv: bool = True
    generate_xlsx: bool = True

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 规范化输入路径|Normalize input path
        self.normalized_input_path = os.path.normpath(os.path.abspath(self.input_path))

        # 发现输入文件|Discover input files
        self.input_files = self._discover_files()

    def _discover_files(self) -> List[str]:
        """
        发现输入文件|Discover input files

        Returns:
            list: 文件路径列表|List of file paths
        """
        path = Path(self.normalized_input_path)

        if path.is_file():
            # 单个文件|Single file
            return [str(path)]
        elif path.is_dir():
            # 目录，查找所有匹配的文件|Directory, find all matching files
            # 支持的文件扩展名|Supported file extensions
            extensions = ['*.fa', '*.fasta', '*.fna', '*.fq', '*.fastq']
            files = []

            for ext in extensions:
                pattern = str(path / ext)
                files.extend(glob.glob(pattern))

            # 同时支持大写扩展名|Also support uppercase extensions
            for ext in ['*.FA', '*.FASTA', '*.FNA', '*.FQ', '*.FASTQ']:
                pattern = str(path / ext)
                files.extend(glob.glob(pattern))

            return sorted(list(set(files)))  # 去重并排序|Remove duplicates and sort
        else:
            raise ValueError(f"输入路径不存在|Input path does not exist: {self.normalized_input_path}")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入路径|Check input path
        if not os.path.exists(self.normalized_input_path):
            errors.append(f"输入路径不存在|Input path does not exist: {self.normalized_input_path}")

        # 检查是否找到文件|Check if files were found
        if not self.input_files:
            errors.append(f"未找到有效的FASTA/FASTQ文件|No valid FASTA/FASTQ files found in: {self.normalized_input_path}")

        # 检查参数范围|Check parameter ranges
        if self.min_length < 0:
            errors.append(f"最小长度必须为非负整数|Minimum length must be non-negative: {self.min_length}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
