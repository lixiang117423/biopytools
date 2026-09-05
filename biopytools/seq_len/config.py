"""seq_len 配置类|seq_len configuration dataclass"""

import os
from dataclasses import dataclass
from typing import List, Optional

from ..common.paths import expand_path

# 视为 FASTA 的扩展名(小写)|extensions treated as FASTA (lowercase)
_FASTA_EXTS = ('.fa', '.fasta', '.fna', '.faa', '.ffn', '.frn')


def _is_fasta_name(name: str) -> bool:
    """文件名是否像 FASTA(支持 .gz)|Does filename look like FASTA (.gz supported)"""
    low = name.lower()
    if low.endswith('.gz'):
        low = low[:-3]
    return low.endswith(_FASTA_EXTS)


def _looks_like_file(path: str) -> bool:
    """路径是否更像文件(含扩展名且非已存在目录)|Heuristic: does path look like a file?"""
    if os.path.isdir(path):
        return False
    return ('.' in os.path.basename(path)) and (not path.endswith(os.sep))


def _summary_path_from_tsv(tsv_path: str) -> str:
    """由主表路径推导汇总表路径(X.tsv → X_summary.tsv)|Derive summary path from main tsv"""
    root, ext = os.path.splitext(tsv_path)
    return f"{root}_summary{ext}"


@dataclass
class SeqLenConfig:
    """FASTA 序列长度统计配置|FASTA sequence length statistics config"""

    input_path: str
    output: str
    prefix: Optional[str] = None
    min_length: int = 0
    sort: bool = False
    summary: bool = True
    log_file: Optional[str] = None
    log_level: str = 'INFO'
    verbose: bool = False

    def __post_init__(self):
        """展开路径 + 解析输入(文件/文件夹)+ 智能 -o|Expand paths, resolve input, smart -o"""
        # 先展开输入再判断文件/文件夹(~ 需先展开)|expand input before isdir check
        self.input_path = expand_path(self.input_path)
        self.is_folder = os.path.isdir(self.input_path)
        if self.log_file:
            self.log_file = expand_path(self.log_file)

        # 解析输入文件列表|resolve input file list
        if self.is_folder:
            self.input_files: List[str] = sorted(
                os.path.join(self.input_path, n)
                for n in os.listdir(self.input_path) if _is_fasta_name(n))
        else:
            self.input_files = [self.input_path]

        self._resolve_output()

    def _derive_prefix(self) -> str:
        """从输入路径推断 prefix(文件→stem,文件夹→目录名)|Derive prefix from input (file stem / dir name)"""
        base = os.path.basename(self.input_path.rstrip(os.sep))
        low = base.lower()
        if low.endswith('.gz'):  # 先剥 .gz|strip .gz first
            base = base[:-3]
            low = low[:-3]
        for ext in _FASTA_EXTS:
            if low.endswith(ext):
                return base[:-len(ext)]
        return base  # 目录名或无 FASTA 扩展名|dir name or no fasta ext

    def _resolve_output(self):
        """智能 -o:目录(写 {prefix}_seq_len.tsv)vs 文件(直接当主表)|Smart -o: dir vs file"""
        raw = self.output
        expanded = expand_path(raw)
        is_dir_target = (raw.endswith(os.sep) or raw.endswith('/')
                         or os.path.isdir(expanded)
                         or not _looks_like_file(raw))

        if is_dir_target:
            self.prefix = self.prefix or self._derive_prefix()
            output_dir = expand_path(raw.rstrip('/').rstrip(os.sep) or '.')
            os.makedirs(output_dir, exist_ok=True)
            self.tsv_path = os.path.join(output_dir, f"{self.prefix}_seq_len.tsv")
        else:
            self.tsv_path = expanded
            output_dir = os.path.dirname(expanded) or '.'
            os.makedirs(output_dir, exist_ok=True)
            # 文件模式 prefix 仅用于日志,默认从输入推导|file-mode prefix for logging only
            self.prefix = self.prefix or self._derive_prefix()

        self.output_dir = output_dir
        self.summary_path = _summary_path_from_tsv(self.tsv_path)

    def validate(self):
        """校验输入与参数|Validate input and parameters"""
        errors = []
        if not os.path.exists(self.input_path):
            errors.append(f"输入路径不存在|Input path not found: {self.input_path}")
        if self.is_folder and not self.input_files:
            errors.append(f"目录下无 FASTA 文件|No FASTA file in directory: {self.input_path}")
        if self.min_length < 0:
            errors.append("min_length 不能为负|must be >= 0")
        if errors:
            raise ValueError("\n".join(errors))
        return True
