"""
JCVI共享配置基类|JCVI Shared Configuration Base Class
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from ..common.paths import expand_path


@dataclass
class JcviBaseConfig:
    """JCVI基础配置类|JCVI Base Configuration

    mcscan和allelic共享的配置字段, 各子命令可扩展此基类
    Shared config fields for mcscan and allelic, subcommands extend this base
    """

    # 必需参数|Required parameters
    input_dir: str = ""
    output_dir: str = ""

    # conda环境|conda environment
    conda_env: str = "JCVI_v.1.5.6"

    # 输入文件扩展名|Input file extensions
    gff_ext: str = ".gff"
    fa_ext: str = ".fa"
    pep_ext: str = ".pep"

    # JCVI参数|JCVI parameters
    gff_type: str = "mRNA"
    gff_key: str = "ID"
    cscore: float = 0.7
    min_size: int = 4
    dbtype: str = "prot"
    align_soft: str = "last"
    no_strip_names: bool = True
    no_dotplot: bool = True

    # 性能参数|Performance parameters
    threads: int = 24

    # 指定比较配对, 格式: ["A,B", "A,C"]; None=全两两|Specific pairs; None=all pairwise
    pairs: Optional[List[str]] = None

    def __post_init__(self):
        self.input_dir = expand_path(self.input_dir)
        self.output_dir = expand_path(self.output_dir)
        os.makedirs(self.output_dir, exist_ok=True)

    def validate(self):
        errors = []
        if not self.input_dir or not Path(self.input_dir).is_dir():
            errors.append(f"输入目录不存在|Input directory not found: {self.input_dir}")
        if not self.output_dir:
            errors.append("输出目录不能为空|Output directory cannot be empty")
        if self.dbtype not in ("prot", "nucl"):
            errors.append(f"dbtype必须为prot或nucl: {self.dbtype}")
        if self.align_soft == "diamond_blastp" and self.dbtype != "prot":
            errors.append("diamond_blastp仅支持prot模式")
        if errors:
            raise ValueError("\n".join(errors))
        return True
