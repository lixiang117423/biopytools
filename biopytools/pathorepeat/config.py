"""pathorepeat配置模块|pathorepeat Configuration Module"""

import os
from dataclasses import dataclass
from typing import Optional

from biopytools.common.paths import expand_path, get_domain_tool_path

# 支持的基因组后缀(未压缩)|Supported genome suffixes (uncompressed)
FASTA_SUFFIXES = ('.fa', '.fna', '.fasta')
# 压缩后缀一律拒绝(明确报错,不做静默解压)|Compressed suffixes are always rejected
FASTA_GZ_SUFFIXES = ('.fa.gz', '.fna.gz', '.fasta.gz')

VALID_LOG_LEVELS = ('DEBUG', 'INFO', 'WARNING', 'ERROR')

# TEsorter -db 可选值,实测 repeat 环境 TEsorter 1.5.1 -h
# |TEsorter -db choices, verified from `TEsorter -h` (repeat env, v1.5.1)
# 注意:REXdb 系以植物/动物 lineage 为主,卵菌/原生生物覆盖弱,可试 gydb
# |Note: REXdb lineages are plant/metazoa-heavy; try gydb for oomycetes/protists
TESORTER_DBS = ('gydb', 'rexdb', 'rexdb-plant', 'rexdb-metazoa', 'rexdb-v3',
                'rexdb-plantv3', 'rexdb-metazoav3', 'rexdb-pnas', 'rexdb-line',
                'sine')

# 屏蔽模式→RepeatMasker 参数(xsmall=小写软屏蔽,help 原文已核实)
# |Masking mode → RepeatMasker flags (xsmall=lowercase soft mask, verified)
# 注意:RM 4.2.4 已移除 -soft,-xsmall 是其现代等价(小写软屏蔽),与 repeatmask
# 模块一致;硬编码 -soft 会让 RepeatMasker 直接报 unknown option
# |Note: RM 4.2.4 dropped -soft; -xsmall is its modern equivalent (lowercase soft
# mask), consistent with the repeatmask module; hardcoding -soft makes RM abort
MASKING_FLAGS = {
    'xsmall': ['-xsmall'],
    # soft 映射为 -xsmall:RM 4.2.4 已移除 -soft 选项(实测 "Unknown option: soft"),
    # 小写软屏蔽语义由 -xsmall 承担(与 repeatmask 模块处理一致)
    # |soft maps to -xsmall: RM 4.2.4 dropped the -soft flag (verified);
    # lowercase soft-masking is carried by -xsmall
    'soft': ['-xsmall'],
    'hard': [],
    'x': ['-x'],
}


@dataclass
class PathorepeatConfig:
    """pathorepeat配置类|pathorepeat Configuration Class

    单 FASTA=单样本;文件夹=批量(逐样品跑完整四步)
    |Single FASTA=one sample; directory=batch (full 4 steps per sample)
    """

    input: str = ''
    output_dir: str = './pathorepeat_output'
    threads: int = 12
    masking_mode: str = 'xsmall'
    ltr_struct: bool = True
    tesorter_db: str = 'rexdb'
    db_hmm: Optional[str] = None
    # Dfam famdb 数据目录(含 famdb.py 与 *.h5);设置后注入 FAMDB_DIR 环境变量,
    # 使 RepeatModeler 分类可用 RM2 自带 Dfam 参考
    # |Dfam famdb data dir (famdb.py + *.h5); injected as FAMDB_DIR so the
    # RepeatModeler classification step can use its built-in Dfam reference
    famdb_dir: Optional[str] = None
    effector_bed: Optional[str] = None
    effector_gff: Optional[str] = None
    genome_name: Optional[str] = None
    skip_completed: bool = True
    log_level: str = 'INFO'

    # 软件路径(env_map 注册 repeat 域环境)|Tool paths (repeat domain via env_map)
    build_database_path: str = ''
    repeatmodeler_path: str = ''
    repeatmasker_path: str = ''
    tesorter_path: str = ''

    def __post_init__(self):
        """初始化后处理:~ 展开 + 域环境路径解析|Post-init: expand ~ + domain paths"""
        for attr in ('input', 'output_dir', 'db_hmm', 'effector_bed', 'effector_gff',
                     'famdb_dir'):
            val = getattr(self, attr)
            if val:
                setattr(self, attr, expand_path(val))
        if not self.build_database_path:
            self.build_database_path = get_domain_tool_path(
                'BuildDatabase', '~/miniforge3/envs/repeat/bin/BuildDatabase',
                'BUILDDATABASE_PATH')
        if not self.repeatmodeler_path:
            self.repeatmodeler_path = get_domain_tool_path(
                'RepeatModeler', '~/miniforge3/envs/repeat/bin/RepeatModeler',
                'REPEATMODELER_PATH')
        if not self.repeatmasker_path:
            self.repeatmasker_path = get_domain_tool_path(
                'RepeatMasker', '~/miniforge3/envs/repeat/bin/RepeatMasker',
                'REPEATMASKER_PATH')
        if not self.tesorter_path:
            self.tesorter_path = get_domain_tool_path(
                'TEsorter', '~/miniforge3/envs/repeat/bin/TEsorter',
                'TESORTER_PATH')
        # 显式传入的 ~ 工具路径同样必须展开,否则 os.path.exists 恒 False
        # |Explicitly-passed ~ tool paths must expand too, else exists() is always False
        for attr in ('build_database_path', 'repeatmodeler_path',
                     'repeatmasker_path', 'tesorter_path'):
            setattr(self, attr, expand_path(getattr(self, attr)))
        os.makedirs(self.output_dir, exist_ok=True)

    @property
    def is_batch(self) -> bool:
        """文件夹输入=批量模式|Directory input means batch mode"""
        return bool(self.input) and os.path.isdir(self.input)

    def sample_name(self, path: str) -> str:
        """样品名:单文件模式尊重 --genome-name|Sample name honors --genome-name"""
        if not self.is_batch and self.genome_name:
            return self.genome_name
        from .utils import genome_name
        return genome_name(path)

    def validate(self):
        """验证配置(收集全部错误一次抛出)|Validate (collect all errors, raise once)"""
        errors = []

        # 输入存在性|Input existence
        if not self.input:
            errors.append("必须提供输入|-i/--input is required")
        elif not os.path.exists(self.input):
            errors.append(f"输入路径不存在|Input path not found: {self.input}")
        elif os.path.isfile(self.input):
            name = os.path.basename(self.input).lower()
            if name.endswith(FASTA_GZ_SUFFIXES):
                errors.append(f"不支持压缩FASTA,请先解压|Compressed FASTA not supported, "
                              f"decompress first: {self.input}")
            elif not name.endswith(FASTA_SUFFIXES):
                errors.append(f"输入应为 {', '.join(FASTA_SUFFIXES)} FASTA"
                              f"|Input must be a FASTA: {self.input}")

        # 数值与选择项|Numerics and choices
        if self.threads < 1:
            errors.append(f"线程数必须为正|Threads must be positive: {self.threads}")
        if self.masking_mode not in MASKING_FLAGS:
            errors.append(f"无效屏蔽模式|Invalid masking mode: {self.masking_mode} "
                          f"(可选|choices: {', '.join(MASKING_FLAGS)})")
        if self.tesorter_db not in TESORTER_DBS:
            errors.append(f"无效TEsorter数据库|Invalid TEsorter db: {self.tesorter_db} "
                          f"(可选|choices: {', '.join(TESORTER_DBS)})")
        if self.log_level.upper() not in VALID_LOG_LEVELS:
            errors.append(f"日志级别无效|Invalid log level: {self.log_level}")

        # effector 互斥 + 批量模式禁用|effector exclusivity + batch restrictions
        if self.effector_bed and self.effector_gff:
            errors.append("--effector-bed 与 --effector-gff 互斥"
                          "|--effector-bed and --effector-gff are exclusive")
        if self.is_batch:
            forbidden = [n for n, v in (('effector_bed', self.effector_bed),
                                        ('effector_gff', self.effector_gff),
                                        ('genome_name', self.genome_name)) if v]
            if forbidden:
                errors.append(f"批量(文件夹)模式不支持参数|Batch (directory) mode does not "
                              f"support: {', '.join('--' + n.replace('_', '-') for n in forbidden)}"
                              f"(effector 交叉检查属单样本分析,多基因组 seqid 会撞名"
                              f"|effector cross-check is single-sample only)")

        # 附属文件存在性|Aux file existence
        for label, path in (('--db-hmm', self.db_hmm),
                            ('--effector-bed', self.effector_bed),
                            ('--effector-gff', self.effector_gff)):
            if path and not os.path.exists(path):
                errors.append(f"{label} 文件不存在|File not found: {path}")
        if self.famdb_dir:
            if not os.path.isdir(self.famdb_dir):
                errors.append(f"--famdb-dir 目录不存在|Directory not found: "
                              f"{self.famdb_dir}")
            elif not os.path.exists(os.path.join(self.famdb_dir, 'famdb.py')):
                errors.append(f"--famdb-dir 下未找到 famdb.py(可将 conda 环境内的 "
                              f"famdb.py 软链到该目录)|famdb.py not found in --famdb-dir "
                              f"(symlink the one from the conda env): {self.famdb_dir}")

        # 工具二进制|Tool binaries
        for attr in ('build_database_path', 'repeatmodeler_path',
                     'repeatmasker_path', 'tesorter_path'):
            if not os.path.exists(getattr(self, attr)):
                errors.append(f"工具未找到|Tool not found: {getattr(self, attr)} "
                              f"(检查repeat环境或对应*_PATH环境变量"
                              f"|check repeat env or the *_PATH env var)")

        if errors:
            raise ValueError("\n".join(errors))
        return True
