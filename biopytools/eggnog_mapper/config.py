"""
eggnog-mapper (emapper) 功能注释配置模块|eggnog-mapper functional annotation config module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path, get_tool_path


VALID_ITYPES = ("proteins", "CDS", "genome", "metagenome")
VALID_MODES = ("mmseqs", "diamond", "hmmer", "no_search", "cache")
VALID_FASTA_EXTS = (".fa", ".faa", ".fasta", ".ffn", ".fna")


@dataclass
class EggnogMapperConfig:
    """eggnog-mapper 配置类|eggnog-mapper configuration class"""

    # 必需参数|Required
    input_file: str
    output_dir: str

    # 可选参数|Optional
    itype: str = "proteins"
    translate: bool = False
    mode: str = "mmseqs"
    cpu: int = 12
    sensmode: str = "sensitive"
    seed_ortholog_evalue: float = 0.001
    prefix: Optional[str] = None
    resume: bool = False
    override: bool = False
    no_format: bool = False

    # 路径(默认走 get_tool_path,支持~与env)|Paths via get_tool_path
    emapper_path: str = field(
        default_factory=lambda: get_tool_path(
            "emapper",
            "~/miniforge3/envs/eggnog-mapper_v.2.1.15/bin/emapper.py",
            "EMAPPER_PATH",
        )
    )
    data_dir: str = field(
        default_factory=lambda: get_tool_path(
            "eggnog_db", "~/database/eggnog", "EGGNOG_DATA_DIR"
        )
    )

    def __post_init__(self):
        """初始化后处理:展开路径、建目录、校验枚举|Expand paths, create dirs, validate enums."""
        # 展开路径(规范§11.3.1)|Expand paths
        self.input_file = expand_path(self.input_file)
        self.output_dir = expand_path(self.output_dir)
        self.emapper_path = expand_path(self.emapper_path)
        self.data_dir = expand_path(self.data_dir)

        # 创建输出目录结构|Create output structure
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        for sub in ("00_pipeline_info", "01_emapper", "99_logs"):
            (self.output_path / sub).mkdir(parents=True, exist_ok=True)

        # 前缀默认=输入 stem|Prefix default = input stem
        self.prefix = self.prefix or Path(self.input_file).stem

        # 枚举校验|Enum validation
        if self.itype not in VALID_ITYPES:
            raise ValueError(
                f"无效的itype|Invalid itype: {self.itype} "
                f"(必须是|must be one of: {', '.join(VALID_ITYPES)})"
            )
        if self.mode not in VALID_MODES:
            raise ValueError(
                f"无效的mode|Invalid mode: {self.mode} "
                f"(必须是|must be one of: {', '.join(VALID_MODES)})"
            )

    def validate(self):
        """校验输入/工具/数据库就绪状态|Validate input/tool/DB readiness."""
        errors = []

        # 输入文件|Input file
        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")
        else:
            ext = os.path.splitext(self.input_file)[1].lower()
            if ext and ext not in VALID_FASTA_EXTS:
                errors.append(
                    f"输入文件应为FASTA格式|Input should be FASTA: {self.input_file}"
                )

        # emapper 可执行|emapper executable
        if not os.path.exists(self.emapper_path):
            errors.append(
                f"emapper.py不存在|emapper.py not found: {self.emapper_path}"
            )

        # 数据库目录及必需库|DB dir and required files
        if not os.path.isdir(self.data_dir):
            errors.append(
                f"数据库目录不存在|DB dir not found: {self.data_dir}\n"
                f"  提示|Hint: 用|use download_eggnog_data.py 下载|download:\n"
                f"    conda run -n eggnog-mapper_v.2.1.15 --no-capture-output "
                f"download_eggnog_data.py -y -M --data_dir {self.data_dir}"
            )
        else:
            for required in ("eggnog.db", "eggnog.taxa.db"):
                if not os.path.exists(os.path.join(self.data_dir, required)):
                    errors.append(
                        f"缺少必需库|Missing required DB: "
                        f"{self.data_dir}/{required}"
                    )
            if self.mode == "mmseqs":
                names = os.listdir(self.data_dir)
                mmseqs_ok = any(
                    n.startswith("eggnog_mmseqs") or n == "mmseqs" for n in names
                )
                if not mmseqs_ok:
                    errors.append(
                        f"缺少mmseqs库|Missing mmseqs DB in {self.data_dir}\n"
                        f"  提示|Hint: download_eggnog_data.py -M --data_dir "
                        f"{self.data_dir}"
                    )
            if self.mode == "diamond":
                dmnd = os.path.join(self.data_dir, "eggnog_proteins.dmnd")
                if not os.path.exists(dmnd):
                    errors.append(
                        f"缺少diamond库|Missing eggnog_proteins.dmnd: {dmnd}\n"
                        f"  提示|Hint: download_eggnog_data.py -y --data_dir "
                        f"{self.data_dir}"
                    )

        if errors:
            raise ValueError("\n".join(errors))
        return True
