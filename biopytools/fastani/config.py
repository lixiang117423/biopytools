"""fastANI配置模块|fastANI Configuration Module"""

import os
from dataclasses import dataclass
from typing import Optional

from biopytools.common.paths import expand_path, get_domain_tool_path

# 支持的基因组后缀|Supported genome suffixes
FASTA_SUFFIXES = ('.fa', '.fna', '.fasta')
FASTA_GZ_SUFFIXES = ('.fa.gz', '.fna.gz', '.fasta.gz')

VALID_LOG_LEVELS = ('DEBUG', 'INFO', 'WARNING', 'ERROR')


@dataclass
class FastaniConfig:
    """fastANI配置类|fastANI Configuration Class

    两种互斥模式|Two exclusive modes:
      A: input(目录/列表/单fasta)→ all-vs-all
      B: query + reference → 定向比较|directional comparison
    """

    # 模式A输入|Mode A input
    input: Optional[str] = None
    # 模式B输入|Mode B inputs
    query: Optional[str] = None
    reference: Optional[str] = None

    output_dir: str = './fastani_output'

    # fastANI可调参数|fastANI tunables
    threads: int = 12
    kmer: int = 16
    frag_len: int = 3000
    min_fraction: float = 0.2

    log_level: str = 'INFO'

    # 软件路径|Software path(env_map注册pop域环境)
    fastani_path: str = ''

    def __post_init__(self):
        """初始化后处理|Post-initialization"""
        for attr in ('input', 'query', 'reference', 'output_dir'):
            val = getattr(self, attr)
            if val:
                setattr(self, attr, expand_path(val))
        if not self.fastani_path:
            self.fastani_path = get_domain_tool_path(
                'fastANI', '~/miniforge3/envs/pop/bin/fastANI', 'FASTANI_PATH')
        os.makedirs(self.output_dir, exist_ok=True)

    @property
    def all_vs_all(self) -> bool:
        """是否all-vs-all模式|Whether all-vs-all mode"""
        return self.input is not None

    def validate(self):
        """验证配置(收集全部错误一次抛出)|Validate (collect all errors, raise once)"""
        errors = []

        # 模式互斥|Mode exclusivity
        has_a = self.input is not None
        has_b = self.query is not None or self.reference is not None
        if has_a and has_b:
            errors.append("输入互斥|Inputs are exclusive: -i 与 -q/-r 不可同时使用|"
                          "-i and -q/-r cannot be combined")
        if not has_a and not has_b:
            errors.append("至少提供一种输入|Provide at least one input: "
                          "-i(all-vs-all)或 -q+-r(定向)|-i (all-vs-all) or -q+-r (directional)")
        if has_b and not (self.query is not None and self.reference is not None):
            errors.append("定向模式必须同时提供|Directional mode requires both: "
                          "-q 与 reference|-q and reference")

        # 路径存在性|Path existence
        for label, path in (("输入|Input", self.input),
                            ("query", self.query), ("reference", self.reference)):
            if path and not os.path.exists(path):
                errors.append(f"{label}路径不存在|Path not found: {path}")

        # fastANI 二进制|Binary
        if not os.path.exists(self.fastani_path):
            errors.append(f"fastANI未找到|fastANI not found: {self.fastani_path} "
                          "(检查pop环境或设置FASTANI_PATH|check pop env or set FASTANI_PATH)")

        # 数值边界|Numeric bounds
        if self.threads < 1:
            errors.append(f"线程数必须为正|Threads must be positive: {self.threads}")
        if not 1 <= self.kmer <= 16:
            errors.append(f"kmer必须在1-16之间|kmer must be 1-16: {self.kmer}")
        if self.frag_len < 1:
            errors.append(f"片段长度必须为正|frag_len must be positive: {self.frag_len}")
        if not 0 < self.min_fraction <= 1:
            errors.append(f"min_fraction必须在(0,1]|min_fraction must be in (0,1]: "
                          f"{self.min_fraction}")
        if self.log_level.upper() not in VALID_LOG_LEVELS:
            errors.append(f"日志级别无效|Invalid log level: {self.log_level}")

        if errors:
            raise ValueError("\n".join(errors))
        return True
