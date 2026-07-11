"""gene_table 配置类|gene_table configuration dataclass"""

import os
from dataclasses import dataclass, field
from typing import Optional, Tuple
from ..common.paths import expand_path, get_tool_path

# 视为输出的 GFF 扩展名(用于从文件名推断 prefix)|GFF extensions stripped when deriving prefix
_GFF_EXT = ('.gff3.gz', '.gff.gz', '.gff3', '.gff')
# 表文件扩展名(用于剥离出 prefix)|table extensions stripped to derive prefix
_TABLE_EXT = ('.gene_table.tsv', '.tsv')
# 无样本语义的通用文件名(回退到 GFF 推断 prefix)|generic names that fall back to GFF-derived prefix
_GENERIC_NAMES = ('gene_table', 'output', 'out', 'result', 'results')


def _derive_prefix_from_gff(gff_file: str) -> str:
    """从 GFF 文件名推断 prefix|Derive prefix from GFF basename"""
    base = os.path.basename(gff_file)
    for ext in _GFF_EXT:
        if base.lower().endswith(ext):
            return base[:-len(ext)]
    return os.path.splitext(base)[0]


def _looks_like_file(path: str) -> bool:
    """路径是否更像文件(含扩展名且非已存在目录)|Heuristic: does path look like a file?"""
    if os.path.isdir(path):
        return False
    return ('.' in os.path.basename(path)) and (not path.endswith(os.sep))


@dataclass
class GeneTableConfig:
    """基因信息+序列合并表配置|Gene info + sequence merged table config"""

    genome_file: str
    gff_file: str
    output: str
    prefix: Optional[str] = None
    longest_only: bool = False
    transcript_types: Tuple[str, ...] = ('mRNA', 'transcript')
    gene_type: str = 'gene'
    min_length: int = 0
    gffread_path: str = field(
        default_factory=lambda: get_tool_path('gffread', '~/.local/bin/gffread', 'GFFREAD_PATH'))
    log_file: Optional[str] = None
    log_level: str = 'INFO'
    verbose: bool = False

    def __post_init__(self):
        """初始化后处理:展开路径 + 智能 -o 解析|Post-init: expand paths + smart -o resolution"""
        # 关键:展开所有含 ~ 的路径|CRITICAL: expand all ~ paths
        self.genome_file = expand_path(self.genome_file)
        self.gff_file = expand_path(self.gff_file)
        self.gffread_path = expand_path(self.gffread_path)
        if self.log_file:
            self.log_file = expand_path(self.log_file)

        prefix_given = self.prefix is not None

        # 智能 -o 解析:目录 vs 文件|Smart -o resolution (directory vs file)
        if self.output.endswith(os.sep) or (os.path.isdir(self.output) and not _looks_like_file(self.output)):
            # 当作输出目录|Treat as output directory
            self.output_dir = self.output.rstrip(os.sep) or '.'
            self.prefix = self.prefix or _derive_prefix_from_gff(self.gff_file)
            self.tsv_path = os.path.join(self.output_dir, f"{self.prefix}.gene_table.tsv")
        else:
            # 当作表文件路径|Treat as table file path
            self.tsv_path = self.output
            self.output_dir = os.path.dirname(self.output) or '.'
            if not prefix_given:
                base = os.path.basename(self.output)
                stripped = base
                for ext in _TABLE_EXT:
                    if stripped.lower().endswith(ext):
                        stripped = stripped[:-len(ext)]
                        break
                # 通用文件名回退到 GFF 推断|Generic names fall back to GFF-derived prefix
                if stripped and stripped.lower() not in _GENERIC_NAMES:
                    self.prefix = stripped
                else:
                    self.prefix = _derive_prefix_from_gff(self.gff_file)
            # prefix_given 时 self.prefix 已是用户值,保持不变|user prefix kept as-is

        os.makedirs(self.output_dir, exist_ok=True)
        self.gene_fa = os.path.join(self.output_dir, f"{self.prefix}.gene.fa")
        self.cds_fa = os.path.join(self.output_dir, f"{self.prefix}.cds.fa")
        self.pep_fa = os.path.join(self.output_dir, f"{self.prefix}.pep.fa")

    def validate(self):
        """校验输入文件存在|Validate that input files exist"""
        errors = []
        for label, path in (('基因组|genome', self.genome_file),
                            ('GFF', self.gff_file),
                            ('gffread', self.gffread_path)):
            if not os.path.exists(path):
                errors.append(f"{label} 不存在|not found: {path}")
        if self.min_length < 0:
            errors.append("min_length 不能为负|must be >= 0")
        if errors:
            raise ValueError("\n".join(errors))
