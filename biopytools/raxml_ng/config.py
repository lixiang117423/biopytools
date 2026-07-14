"""
RAxML-NG 系统发育树分析配置模块|RAxML-NG Phylogenetic Tree Analysis Configuration Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import get_tool_path, expand_path


@dataclass
class RaxmlNgConfig:
    """RAxML-NG 系统发育树分析配置类|RAxML-NG Phylogenetic Tree Analysis Configuration Class"""

    # 必需参数|Required parameters
    input_file: str
    output_dir: str

    # 核心参数|Core parameters
    prefix: Optional[str] = None  # None 表示用输入文件名|None means use input filename
    mode: str = 'all'  # all / search / support
    model: Optional[str] = None  # None 表示自适应选择|None means adaptive selection
    threads: int = 12

    # Bootstrap 参数|Bootstrap parameters
    bs_trees: str = '1000'  # 重复数或 autoMRE{N};all/support 用|replicates or autoMRE{N}
    bs_metric: str = 'fbp'  # fbp/tbe/sh/ebg/rbs/ps/pbs/ic1/ica/gcf

    # 树选项|Tree options
    tree: Optional[str] = None  # search 起始树 / support 参考树|search starting / support reference
    outgroup: Optional[str] = None

    # 随机与续传|Random and resume
    seed: Optional[int] = None
    redo: bool = False

    # 工具路径|Tool path (静态二进制,不在 conda /envs/ 下,直接调用)|Static binary, not under conda /envs/
    raxml_ng_path: str = field(
        default_factory=lambda: get_tool_path(
            'raxml_ng',
            '~/software/RAxML/raxml-ng_v2.0.2_linux_x86_64/raxml-ng',
            'RAXML_NG_PATH'
        )
    )

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开工具路径(关键:必须展开~)|Expand tool path (CRITICAL: must expand ~)
        self.raxml_ng_path = expand_path(self.raxml_ng_path)

        # 输出目录及规范子目录|Output dir and standard subdirs
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.logs_path = self.output_path / '99_logs'
        self.logs_path.mkdir(parents=True, exist_ok=True)
        self.info_path = self.output_path / '00_pipeline_info'
        self.info_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.input_file = os.path.normpath(os.path.abspath(self.input_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 前缀默认用输入文件名(去扩展名)|Prefix defaults to input filename stem
        if self.prefix is None:
            self.prefix = Path(self.input_file).stem

        if self.tree:
            self.tree = os.path.normpath(os.path.abspath(self.tree))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file does not exist: {self.input_file}")

        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread number must be positive: {self.threads}")

        if self.mode not in ('all', 'search', 'support'):
            errors.append(f"模式无效|Invalid mode (all/search/support): {self.mode}")

        valid_metrics = ('fbp', 'tbe', 'sh', 'ebg', 'rbs', 'ps', 'pbs', 'ic1', 'ica', 'gcf')
        if self.bs_metric not in valid_metrics:
            errors.append(f"支持值类型无效|Invalid bs-metric {self.bs_metric}: {valid_metrics}")

        # support 模式需要参考树和 bootstrap 树文件|support mode requires reference tree + bootstrap tree file
        if self.mode == 'support':
            if not self.tree:
                errors.append("support 模式需要参考树(--tree)|support mode requires reference tree (--tree)")
            elif not os.path.exists(self.tree):
                errors.append(f"参考树文件不存在|Reference tree file does not exist: {self.tree}")
            if not os.path.exists(self.bs_trees):
                errors.append(
                    f"support 模式 --bs-trees 需为 bootstrap 树文件|"
                    f"support mode --bs-trees must be a bootstrap tree file: {self.bs_trees}"
                )

        # search 模式的起始树(可选)|search starting tree (optional)
        if self.tree and self.mode != 'support' and not os.path.exists(self.tree):
            errors.append(f"起始树文件不存在|Starting tree file does not exist: {self.tree}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
