"""
FAPROTAX配置管理模块|FAPROTAX Configuration Management Module
"""

import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path, get_tool_path


VALID_NORMALIZE = {
    'none', 'columns_before_collapsing', 'rows_before_collapsing',
    'columns_after_collapsing', 'rows_after_collapsing',
    'columns_before_collapsing_excluding_unassigned',
    'rows_before_collapsing_excluding_unassigned'
}

VALID_OUTPUT_FORMATS = {'auto', 'BIOM', 'classical'}

VALID_AVERAGE = {
    'none', 'across_records', 'across_group_members',
    'across_used_group_members', 'maximum', 'minimum',
    'minimum_across_records'
}


@dataclass
class FaprotaxtaxConfig:
    """FAPROTAX配置类|FAPROTAX Configuration Class"""

    # ===== 必需参数|Required parameters =====
    input_table: str

    # ===== 输出配置|Output configuration =====
    output_dir: str = "./faprotaxtax_output"

    # ===== 工具路径配置|Tool path configuration =====
    groups_file: str = field(
        default_factory=lambda: get_tool_path(
            'faprotaxtax_groups',
            '~/software/FAPROTAX/FAPROTAX_1.2.12/FAPROTAX.txt',
            'FAPROTAX_GROUPS_PATH'
        )
    )
    collapse_table_path: str = field(
        default_factory=lambda: get_tool_path(
            'collapse_table',
            '~/software/FAPROTAX/FAPROTAX_1.2.12/collapse_table.py',
            'COLLAPSE_TABLE_PATH'
        )
    )
    python_interpreter: str = sys.executable

    # ===== 核心参数|Core parameters =====
    collapse_by_metadata: Optional[str] = None
    group_leftovers_as: Optional[str] = None
    normalize: str = "none"
    average: str = "none"
    row_names_are_in_column: Optional[str] = None
    output_format: str = "auto"

    # ===== 控制参数|Control parameters =====
    threads: int = 1
    force: bool = False
    verbose: bool = False

    # ===== 内部路径（__post_init__中填充）|Internal paths (populated in __post_init__) =====
    output_path: Optional[Path] = field(default=None, init=False, repr=False)
    info_dir: Optional[str] = field(default=None, init=False, repr=False)
    collapsed_dir: Optional[str] = field(default=None, init=False, repr=False)
    report_dir: Optional[str] = field(default=None, init=False, repr=False)
    log_dir: Optional[str] = field(default=None, init=False, repr=False)
    work_dir: Optional[str] = field(default=None, init=False, repr=False)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有路径|Expand all paths
        self.input_table = expand_path(self.input_table)
        self.output_dir = expand_path(self.output_dir)
        self.groups_file = expand_path(self.groups_file)
        self.collapse_table_path = expand_path(self.collapse_table_path)
        self.python_interpreter = expand_path(self.python_interpreter)

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准子目录|Standard subdirectories
        self.info_dir = os.path.join(self.output_dir, "00_pipeline_info")
        self.collapsed_dir = os.path.join(self.output_dir, "01_collapsed")
        self.report_dir = os.path.join(self.output_dir, "02_report")
        self.log_dir = os.path.join(self.output_dir, "99_logs")
        self.work_dir = os.path.join(self.output_dir, "work")

        for dir_path in [self.info_dir, self.collapsed_dir,
                         self.report_dir, self.log_dir, self.work_dir]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.input_table):
            errors.append(f"输入文件不存在|Input file not found: {self.input_table}")

        # 检查组文件|Check groups file
        if not os.path.exists(self.groups_file):
            errors.append(f"功能组数据库不存在|Functional groups database not found: {self.groups_file}")

        # 检查collapse_table.py脚本|Check collapse_table.py script
        if not os.path.exists(self.collapse_table_path):
            errors.append(f"collapse_table.py脚本不存在|collapse_table.py script not found: {self.collapse_table_path}")

        # 验证参数值|Validate parameter values
        if self.normalize not in VALID_NORMALIZE:
            errors.append(f"无效的标准化方式|Invalid normalize: {self.normalize} (可选|options: {VALID_NORMALIZE})")

        if self.output_format not in VALID_OUTPUT_FORMATS:
            errors.append(f"无效的输出格式|Invalid output format: {self.output_format} (可选|options: {VALID_OUTPUT_FORMATS})")

        if self.average not in VALID_AVERAGE:
            errors.append(f"无效的平均方式|Invalid average: {self.average} (可选|options: {VALID_AVERAGE})")

        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
