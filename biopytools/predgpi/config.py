"""
PredGPI配置模块|PredGPI configuration module
"""

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path


# 已知FASTA扩展名(长优先,小写匹配)|Known FASTA extensions (longest first, lowercase match)
_KNOWN_FASTA_EXTS = [
    ".fasta.gz", ".faa.gz", ".fa.gz", ".fna.gz", ".ffn.gz",
    ".fasta", ".faa", ".fa", ".fna", ".ffn",
]


def _derive_prefix(input_file: str) -> str:
    """从输入文件名推导prefix(去已知FASTA扩展)|Derive prefix (strip FASTA ext)"""
    base = os.path.basename(input_file)
    low = base.lower()
    for ext in _KNOWN_FASTA_EXTS:
        if low.endswith(ext):
            return base[:-len(ext)]
    return os.path.splitext(base)[0]


def _sanitize_prefix(prefix: str) -> str:
    """文件名安全化(空格/路径分隔符->_)|Sanitize prefix for filename safety"""
    return re.sub(r"[\s/\\]", "_", prefix)


@dataclass
class PredGPIConfig:
    """PredGPI配置类|PredGPI configuration class"""

    input_file: str
    output_dir: str
    predgpi_home: str = "~/software/predgpi"
    conservative: bool = False
    output_prefix: Optional[str] = None

    def __post_init__(self):
        """初始化后处理:建目录/展路径/推prefix|Post-init: dirs, expand, prefix"""
        # 创建输出目录与99_logs/00_pipeline_info|Create output dir, 99_logs, 00_pipeline_info
        self.output_path = Path(expand_path(self.output_dir))
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.logs_path = self.output_path / "99_logs"
        self.logs_path.mkdir(parents=True, exist_ok=True)
        self.pipeline_info_path = self.output_path / "00_pipeline_info"
        self.pipeline_info_path.mkdir(parents=True, exist_ok=True)

        # 路径标准化与~展开|Normalize and expand ~
        self.input_file = os.path.normpath(os.path.abspath(expand_path(self.input_file)))
        self.output_dir = os.path.normpath(os.path.abspath(expand_path(self.output_dir)))
        self.predgpi_home = os.path.normpath(expand_path(self.predgpi_home))

        # prefix推导(显式优先,否则由文件名)|Derive prefix (explicit wins, else filename)
        raw_prefix = self.output_prefix if self.output_prefix else _derive_prefix(self.input_file)
        self.output_prefix = _sanitize_prefix(raw_prefix)

    def validate(self):
        """验证配置|Validate configuration"""
        errors = []

        if not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")

        if not os.path.isdir(self.predgpi_home):
            errors.append(f"predgpi目录不存在|predgpi home not found: {self.predgpi_home}")

        # 检查必需的模型文件|Check required model files
        gpidat_dir = os.path.join(self.predgpi_home, "GPIDAT")
        if not os.path.isdir(gpidat_dir):
            errors.append(f"GPIDAT目录不存在|GPIDAT directory not found: {gpidat_dir}")
        else:
            # 校验当前模式所需的HMM模型文件+SVM模型|Check HMM+SVM model files for current mode
            hmm_file = "PHMM.TOT.ss.mod_CSDGN" if self.conservative else "PHMM.TOT.ss.mod"
            for mf in (hmm_file, "MOD"):
                if not os.path.isfile(os.path.join(gpidat_dir, mf)):
                    errors.append(f"模型文件不存在|Model file not found: {mf}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
