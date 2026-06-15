"""
RxLR扫描配置管理模块|RxLR Scanner Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class RxLRConfig:
    """RxLR扫描配置类|RxLR Scanner Configuration Class"""

    # 必需参数|Required parameters
    input_file: str
    output_prefix: str

    # 搜索窗口参数|Search window parameters
    window_start: int = 20  # 起始位置(20对应第21位，Python索引|corresponds to position 21, Python index)
    window_end: int = 120   # 结束位置(120对应第120位|corresponds to position 120)
    min_length: int = 120   # 最小序列长度|Minimum sequence length

    # 输出配置|Output configuration
    output_dir: str = "./rxlr_scanner_output"
    excel_output: bool = True
    tsv_output: bool = True

    # 日志配置|Logging configuration
    verbose: bool = False
    log_file: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 生成输出文件名|Generate output filenames
        base_name = self.output_prefix
        self.excel_file = os.path.join(self.output_dir, f"{base_name}.xlsx") if self.excel_output else None
        self.tsv_file = os.path.join(self.output_dir, f"{base_name}.tsv") if self.tsv_output else None
        self.log_file_path = os.path.join(self.output_dir, f"{base_name}.log") if self.log_file else None

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file does not exist: {self.input_file}")

        if not os.path.isfile(self.input_file):
            errors.append(f"输入路径不是文件|Input path is not a file: {self.input_file}")

        # 检查窗口参数|Check window parameters
        if self.window_start < 0:
            errors.append(f"窗口起始位置不能为负数|Window start cannot be negative: {self.window_start}")

        if self.window_end <= self.window_start:
            errors.append(f"窗口结束位置必须大于起始位置|Window end must be greater than start: {self.window_end} <= {self.window_start}")

        if self.min_length < 1:
            errors.append(f"最小序列长度必须大于0|Minimum length must be greater than 0: {self.min_length}")

        # 检查输出配置|Check output configuration
        if not self.excel_output and not self.tsv_output:
            errors.append("至少需要一种输出格式|At least one output format is required (Excel or TSV)")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_output_info(self) -> dict:
        """获取输出文件信息|Get output file information"""
        info = {
            "output_dir": self.output_dir,
        }

        if self.excel_output:
            info["excel_file"] = self.excel_file

        if self.tsv_output:
            info["tsv_file"] = self.tsv_file

        if self.log_file_path:
            info["log_file"] = self.log_file_path

        return info
