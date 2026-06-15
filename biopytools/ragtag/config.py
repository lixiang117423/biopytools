"""
RagTag配置管理模块|RagTag Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class RagTagConfig:
    """RagTag配置类|RagTag Configuration Class"""

    # 必需参数|Required parameters
    reference: str
    query: str
    sample_name: str

    # 可选参数|Optional parameters
    threads: int = 12
    output_dir: str = './ragtag_output'
    prefix: Optional[str] = None  # 序列ID前缀|Sequence ID prefix

    # RagTag scaffolding选项|RagTag scaffolding options
    aligner: str = 'minimap2'  # minimap2, unimap, or nucmer
    min_unique_length: int = 1000  # -f
    min_mapq: int = 10  # -q
    max_merge_distance: int = 100000  # -d
    min_grouping_confidence: float = 0.2  # -i
    min_location_confidence: float = 0.0  # -a
    min_orientation_confidence: float = 0.0  # -s
    concatenate_unplaced: bool = False  # -C
    infer_gaps: bool = False  # -r
    min_gap_size: int = 100  # -g
    max_gap_size: int = 100000  # -m
    remove_small_alignments: bool = False  # --remove-small

    # RagTag软件路径|RagTag software path
    ragtag_path: str = 'ragtag.py'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.reference = os.path.normpath(os.path.abspath(self.reference))
        self.query = os.path.normpath(os.path.abspath(self.query))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 验证aligner参数|Validate aligner parameter
        valid_aligners = ['minimap2', 'unimap', 'nucmer']
        if self.aligner not in valid_aligners:
            raise ValueError(f"不支持的比对器|Unsupported aligner: {self.aligner}. "
                           f"必须是|Must be one of {valid_aligners}")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.reference):
            errors.append(f"参考基因组文件不存在|Reference genome file does not exist: {self.reference}")

        if not os.path.exists(self.query):
            errors.append(f"查询基因组文件不存在|Query genome file does not exist: {self.query}")

        # 检查线程数|Check threads
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        # 检查数值参数范围|Check numeric parameter ranges
        if self.min_unique_length < 0:
            errors.append(f"最小唯一比对长度不能为负数|Min unique alignment length cannot be negative")

        if self.min_mapq < 0:
            errors.append(f"最小MAPQ值不能为负数|Min MAPQ value cannot be negative")

        if not (0.0 <= self.min_grouping_confidence <= 1.0):
            errors.append(f"分组置信度必须在0-1之间|Grouping confidence must be between 0 and 1")

        if not (0.0 <= self.min_location_confidence <= 1.0):
            errors.append(f"位置置信度必须在0-1之间|Location confidence must be between 0 and 1")

        if not (0.0 <= self.min_orientation_confidence <= 1.0):
            errors.append(f"方向置信度必须在0-1之间|Orientation confidence must be between 0 and 1")

        if self.min_gap_size < 0:
            errors.append(f"最小gap大小不能为负数|Min gap size cannot be negative")

        if self.max_gap_size < self.min_gap_size:
            errors.append(f"最大gap大小必须大于等于最小gap大小|Max gap size must be >= min gap size")

        if errors:
            raise ValueError("\n".join(errors))

        return True
