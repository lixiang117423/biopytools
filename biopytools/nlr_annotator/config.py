"""
NLR-Annotator配置管理模块|NLR-Annotator Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path

from ..common.paths import expand_path, get_tool_path


@dataclass
class NLRAnnotatorConfig:
    """NLR-Annotator配置类|NLR-Annotator Configuration Class"""

    # 必需参数|Required parameters
    input_path: str = ""
    output_dir: str = "./output"

    # 目录批处理参数|Directory batch parameters
    sample_suffix: str = "*.cds.fa"

    # 工具路径（通过get_tool_path按优先级获取）|Tool paths (obtained via get_tool_path by priority)
    jar_path: str = ""
    mot_file: str = ""
    store_file: str = ""

    # 运行参数|Runtime parameters
    threads: int = 12
    num_seqs_per_thread: int = 1000

    # 可选输出|Optional outputs
    output_gff: str = ""
    output_bed: str = ""
    output_motifs: str = ""
    output_alignment: str = ""

    # 距离参数|Distance parameters
    distance_within_motif_combination: int = 500
    distance_for_elongating: int = 2500
    distance_between_motif_combinations: int = 50000

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 工具路径：优先使用已设置的值，否则通过get_tool_path获取|Tool paths: use set value first, else get_tool_path
        if not self.jar_path:
            self.jar_path = get_tool_path('nlr_annotator',
                                          '~/software/NLR-Annotator/NLR-Annotator-v2.1b.jar',
                                          'NLR_ANNOTATOR_PATH')
        if not self.mot_file:
            self.mot_file = get_tool_path('nlr_annotator_mot',
                                          '~/software/NLR-Annotator/mot.txt')
        if not self.store_file:
            self.store_file = get_tool_path('nlr_annotator_store',
                                            '~/software/NLR-Annotator/store.txt')

        # 展开所有路径|Expand all paths
        self.jar_path = expand_path(self.jar_path)
        self.mot_file = expand_path(self.mot_file)
        self.store_file = expand_path(self.store_file)
        self.input_path = os.path.abspath(self.input_path)

        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []
        if not self.input_path or not os.path.exists(self.input_path):
            errors.append(f"输入路径不存在|Input path not found: {self.input_path}")
        if not os.path.isfile(self.input_path) and not os.path.isdir(self.input_path):
            errors.append(f"输入路径必须是文件或目录|Input path must be a file or directory: {self.input_path}")
        if not self.jar_path or not os.path.exists(self.jar_path):
            errors.append(f"JAR文件不存在|JAR file not found: {self.jar_path}")
        if not self.mot_file or not os.path.exists(self.mot_file):
            errors.append(f"mot.txt不存在|mot.txt not found: {self.mot_file}")
        if not self.store_file or not os.path.exists(self.store_file):
            errors.append(f"store.txt不存在|store.txt not found: {self.store_file}")
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive")
        if errors:
            raise ValueError("\n".join(errors))
        return True
