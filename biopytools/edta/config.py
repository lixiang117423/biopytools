"""
EDTA转座子注释配置模块|EDTA TE Annotation Configuration Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path, get_domain_tool_path

# EDTA/RepeatMasker 对序列 ID 的长度硬限(超限会导致运行失败)
# |EDTA/RepeatMasker hard limit on sequence-id length (run dies beyond it)
MAX_SEQ_ID_LENGTH = 13

# 默认 EDTA 安装(edta_v.2.3.0 为当前保留独立环境,bin/ 下有 EDTA.pl 与 panEDTA.sh)
# |Default EDTA install (edta_v.2.3.0 is the current legacy standalone env)
DEFAULT_EDTA_ENV = 'edta_v.2.3.0'
DEFAULT_EDTA_PL = f'~/miniforge3/envs/{DEFAULT_EDTA_ENV}/bin/EDTA.pl'
DEFAULT_PANEDTA_SH = f'~/miniforge3/envs/{DEFAULT_EDTA_ENV}/bin/panEDTA.sh'


def _resolve_tool(pl_name: str, default_path: str, explicit_dir: Optional[str],
                  env_var: str) -> Optional[str]:
    """解析 EDTA 工具脚本完整路径|Resolve full path of an EDTA tool script

    优先级|Priority:
      1. 显式安装目录(用户 --edta-path,展开后)|explicit install dir
      2. 功能域/环境解析(get_domain_tool_path:EDTA_PATH 等 env var)|domain resolution
      3. 当前 CONDA_PREFIX/share/EDTA(在 EDTA 环境内运行时)|current CONDA_PREFIX
    """
    if explicit_dir:
        candidate = os.path.join(expand_path(explicit_dir), pl_name)
        if os.path.exists(candidate):
            return candidate
    resolved = get_domain_tool_path(pl_name, default_path, env_var)
    if resolved and os.path.exists(resolved):
        return resolved
    prefix = os.environ.get('CONDA_PREFIX', '')
    if prefix:
        candidate = os.path.join(prefix, 'share', 'EDTA', pl_name)
        if os.path.exists(candidate):
            return candidate
    return None


def _check_sequence_ids(genome: Optional[str], errors: list) -> None:
    """序列 ID 长度预检(防 9 小时级运行中途失败)|Sequence-id length precheck

    EDTA/RepeatMasker 对序列 ID 有 ≤13 字符硬限,超限会在运行数小时后 die;
    此处提前拦截并提示改名工具。
    |EDTA/RepeatMasker dies on ids longer than 13 chars after hours of run;
    intercept early and point to renaming tools.
    """
    if not genome or not os.path.isfile(genome) or genome.lower().endswith('.gz'):
        return
    long_ids = []
    with open(genome, encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            fields = line[1:].split()
            seq_id = fields[0] if fields else ''
            if len(seq_id) > MAX_SEQ_ID_LENGTH:
                long_ids.append(seq_id)
                if len(long_ids) >= 5:
                    break
    if long_ids:
        errors.append(
            f"存在超过{MAX_SEQ_ID_LENGTH}字符的序列ID(EDTA/RepeatMasker 硬限,"
            f"会在运行中失败),示例: {', '.join(long_ids)};"
            f"请先用 biopytools chr_rename 或 assembly_qc 改名后再运行"
            f"|Sequence ids longer than {MAX_SEQ_ID_LENGTH} chars found "
            f"(EDTA/RepeatMasker hard limit), e.g. {', '.join(long_ids)}; "
            f"rename with biopytools chr_rename or assembly_qc first")


@dataclass
class EDTAConfig:
    """EDTA配置类|EDTA Configuration Class"""

    # 必需参数|Required parameters
    genome: str

    # 可选参数|Optional parameters
    species: str = "others"
    step: str = "all"
    overwrite: int = 0
    cds: Optional[str] = None
    curatedlib: Optional[str] = None
    rmlib: Optional[str] = None
    sensitive: int = 0
    anno: int = 0
    rmout: Optional[str] = None
    maxdiv: int = 40
    evaluate: int = 0
    exclude: Optional[str] = None
    force: int = 0
    u: float = 1.3e-8
    threads: int = 12
    debug: int = 0
    output_dir: str = "./edta_output"

    # EDTA路径配置(显式安装目录,可选)|EDTA install dir (explicit, optional)
    edta_path: Optional[str] = None
    repeatmasker_path: Optional[str] = None
    repeatmodeler_path: Optional[str] = None
    ltrretriever_path: Optional[str] = None
    annosine_path: Optional[str] = None

    def __post_init__(self):
        """初始化后处理:~ 展开(§11.B 必须)|Post-init: expand tilde paths"""
        # 全部用户路径先展开再绝对化|Expand tilde before abspath for all inputs
        if self.genome:
            self.genome = os.path.abspath(expand_path(self.genome))
        for attr in ('cds', 'curatedlib', 'rmlib', 'rmout', 'exclude'):
            val = getattr(self, attr)
            if val:
                setattr(self, attr, os.path.abspath(expand_path(val)))
        self.edta_path = expand_path(self.edta_path) if self.edta_path else None
        self.output_dir = expand_path(self.output_dir)

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # step 归一化为小写(API 直调可能传大写)|Normalize step to lowercase
        self.step = str(self.step).lower()

    def _resolve_edta_pl(self) -> Optional[str]:
        """解析 EDTA.pl 完整路径|Resolve full EDTA.pl path"""
        return _resolve_tool('EDTA.pl', DEFAULT_EDTA_PL, self.edta_path,
                             'EDTA_PATH')

    def validate(self):
        """验证配置参数(收集全部错误一次抛出)|Validate (collect all, raise once)"""
        errors = []

        # 必需参数与文件存在性|Required params and file existence
        if not self.genome:
            errors.append("基因组文件未指定|Genome file not specified")
        elif not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome}")

        for label, path in (("CDS文件|CDS", self.cds),
                            ("筛选库文件|Curatedlib", self.curatedlib),
                            ("RepeatModeler库文件|Rmlib", self.rmlib),
                            ("RepeatMasker输出|RMout", self.rmout),
                            ("排除区域文件|Exclude", self.exclude)):
            if path and not os.path.exists(path):
                errors.append(f"{label}不存在|not found: {path}")

        # EDTA.pl 解析(默认走 edta_v.2.3.0 域解析)|Resolve EDTA.pl
        self.resolved_edta_pl = self._resolve_edta_pl()
        if not self.resolved_edta_pl:
            errors.append(
                "未找到 EDTA.pl(检查 edta_v.2.3.0 环境,或设置 --edta-path / "
                "EDTA_PATH 环境变量)|EDTA.pl not found (check the edta_v.2.3.0 "
                "env, or set --edta-path / the EDTA_PATH env var)")

        # 序列 ID 长度预检|Sequence-id length precheck
        _check_sequence_ids(self.genome, errors)

        # 参数范围|Parameter ranges
        if str(self.step).lower() not in ['all', 'filter', 'final', 'anno']:
            errors.append(f"无效的步骤参数|Invalid step parameter: {self.step}")
        if self.overwrite not in [0, 1]:
            errors.append(f"overwrite参数必须为0或1|overwrite must be 0 or 1: {self.overwrite}")
        if self.sensitive not in [0, 1]:
            errors.append(f"sensitive参数必须为0或1|sensitive must be 0 or 1: {self.sensitive}")
        if self.anno not in [0, 1]:
            errors.append(f"anno参数必须为0或1|anno must be 0 or 1: {self.anno}")
        if self.evaluate not in [0, 1]:
            errors.append(f"evaluate参数必须为0或1|evaluate must be 0 or 1: {self.evaluate}")
        if self.force not in [0, 1]:
            errors.append(f"force参数必须为0或1|force must be 0 or 1: {self.force}")
        if self.maxdiv < 0 or self.maxdiv > 100:
            errors.append(f"maxdiv参数必须在0-100之间|maxdiv must be between 0 and 100: {self.maxdiv}")
        if self.threads <= 0:
            errors.append(f"线程数必须大于0|Thread count must be positive: {self.threads}")
        if self.u <= 0:
            errors.append(f"突变率必须大于0|Mutation rate must be positive: {self.u}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class PanEDTAConfig:
    """PanEDTA配置类|PanEDTA Configuration Class"""

    # 必需参数|Required parameters
    genome_list: str

    # 可选参数|Optional parameters
    cds: Optional[str] = None
    curatedlib: Optional[str] = None
    fl_copy: int = 3
    anno: int = 1
    overwrite: int = 0
    threads: int = 12
    output_dir: str = "./panedta_output"

    # EDTA路径配置(显式安装目录,可选)|EDTA install dir (explicit, optional)
    edta_path: Optional[str] = None

    def __post_init__(self):
        """初始化后处理:~ 展开(§11.B 必须)|Post-init: expand tilde paths"""
        if self.genome_list:
            self.genome_list = os.path.abspath(expand_path(self.genome_list))
        if self.cds:
            self.cds = os.path.abspath(expand_path(self.cds))
        if self.curatedlib:
            self.curatedlib = os.path.abspath(expand_path(self.curatedlib))
        self.edta_path = expand_path(self.edta_path) if self.edta_path else None
        self.output_dir = expand_path(self.output_dir)

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def _resolve_panedta_sh(self) -> Optional[str]:
        """解析 panEDTA.sh 完整路径|Resolve full panEDTA.sh path"""
        return _resolve_tool('panEDTA.sh', DEFAULT_PANEDTA_SH, self.edta_path,
                             'PANEDTA_PATH')

    def validate(self):
        """验证配置参数(收集全部错误一次抛出)|Validate (collect all, raise once)"""
        errors = []

        if not self.genome_list:
            errors.append("基因组列表文件未指定|Genome list file not specified")
        elif not os.path.exists(self.genome_list):
            errors.append(f"基因组列表文件不存在|Genome list file not found: {self.genome_list}")

        for label, path in (("CDS文件|CDS", self.cds),
                            ("筛选库文件|Curatedlib", self.curatedlib)):
            if path and not os.path.exists(path):
                errors.append(f"{label}不存在|not found: {path}")

        self.resolved_panedta_sh = self._resolve_panedta_sh()
        if not self.resolved_panedta_sh:
            errors.append(
                "未找到 panEDTA.sh(检查 edta_v.2.3.0 环境,或设置 --edta-path / "
                "PANEDTA_PATH 环境变量)|panEDTA.sh not found (check the "
                "edta_v.2.3.0 env, or set --edta-path / the PANEDTA_PATH env var)")

        if self.anno not in [0, 1]:
            errors.append(f"anno参数必须为0或1|anno must be 0 or 1: {self.anno}")
        if self.overwrite not in [0, 1]:
            errors.append(f"overwrite参数必须为0或1|overwrite must be 0 or 1: {self.overwrite}")
        if self.fl_copy < 0:
            errors.append(f"fl_copy参数必须大于等于0|fl_copy must be >= 0: {self.fl_copy}")
        if self.threads <= 0:
            errors.append(f"线程数必须大于0|Thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
