"""Phytophthora效应子鉴定配置管理模块|Phytophthora Effector Identification Configuration Module"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import get_tool_path, expand_path

_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
DEFAULT_RXLR_HMM = os.path.join(_DATA_DIR, 'paper_RxLR.hmm')
DEFAULT_RXLR_WY_HMM = os.path.join(_DATA_DIR, 'PF18634.hmm')
DEFAULT_CRN_HMM = os.path.join(_DATA_DIR, 'paper_CRN.hmm')
DEFAULT_NLP_HMM = os.path.join(_DATA_DIR, 'paper_nlp.hmm')
DEFAULT_PROTEASE_HMM = os.path.join(_DATA_DIR, 'paper_protease.hmm')
DEFAULT_SCP_HMM = os.path.join(_DATA_DIR, 'paper_scp.hmm')
DEFAULT_ELICITIN_HMM = os.path.join(_DATA_DIR, 'paper_elicitin.hmm')
DEFAULT_YXSL_HMM = os.path.join(_DATA_DIR, 'paper_yxsl.hmm')
DEFAULT_RXLR_BLASTP_QUERIES = os.path.join(_DATA_DIR, 'reference_RxLR_queries.faa')


@dataclass
class PhytoEffectorConfig:
    """Phytophthora效应子鉴定配置类|Phytophthora Effector Identification Configuration Class"""

    mode: str = 'rxlr'
    input_path: str = ''
    output_dir: str = './phyto_effector_output'

    # SignalP参数|SignalP parameters
    skip_signalp: bool = False
    signalp_path: str = field(
        default_factory=lambda: get_tool_path(
            'signalp6', '~/miniforge3/envs/signalp6/bin/signalp6', 'SIGNALP_PATH'
        )
    )
    organism: str = 'eukarya'
    signalp_mode: str = 'slow-sequential'  # fast / slow / slow-sequential
    signalp_version: str = 'both'  # '3', '6', or 'both'
    signalp3_path: str = field(
        default_factory=lambda: get_tool_path(
            'signalp3', '~/miniforge3/envs/signalp_v.3.0b/bin/signalp', 'SIGNALP3_PATH'
        )
    )
    signalp3_sprob_threshold: float = 0.9
    # HMMER参数|HMMER parameters
    hmmsearch_path: str = field(
        default_factory=lambda: get_tool_path(
            'hmmsearch', '~/miniforge3/envs/resistify_v.1.3.0/bin/hmmsearch', 'HMMSEARCH_PATH'
        )
    )

    # BLASTP参数|BLASTP parameters
    blastp_path: str = field(
        default_factory=lambda: get_tool_path(
            'blastp', '~/miniforge3/envs/Blast_v.2.16.0/bin/blastp', 'BLASTP_PATH'
        )
    )
    rxlr_blastp_queries: str = field(
        default_factory=lambda: DEFAULT_RXLR_BLASTP_QUERIES
    )
    rxlr_blastp_evalue: float = 1e-5

    # TMHMM参数|TMHMM parameters
    tmhmm_path: str = field(
        default_factory=lambda: get_tool_path(
            'tmhmm', '~/miniforge3/envs/tmmhmm_v.2.0c/bin/tmhmm', 'TMHMM_PATH'
        )
    )

    # RxLR专属参数|RxLR-specific parameters
    rxlr_hmm: str = DEFAULT_RXLR_HMM
    use_wy_domain: bool = False
    rxlr_wy_hmm: str = DEFAULT_RXLR_WY_HMM
    evalue: float = 1e-5
    score_threshold: float = 0.0

    # CRN专属参数|CRN-specific parameters
    crn_hmm: str = DEFAULT_CRN_HMM

    # 其他效应子类型参数|Other effector type parameters
    nlp_hmm: str = DEFAULT_NLP_HMM
    protease_hmm: str = DEFAULT_PROTEASE_HMM
    scp_hmm: str = DEFAULT_SCP_HMM
    elicitin_hmm: str = DEFAULT_ELICITIN_HMM
    yxsl_hmm: str = DEFAULT_YXSL_HMM

    # 通用参数|Common parameters
    threads: int = 12

    # 内部属性|Internal attributes (不通过__init__传入)
    _input_files: list = field(default_factory=list, init=False)
    _combined_fasta: str = field(default='', init=False)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 展开所有路径|Expand all paths
        self.input_path = expand_path(self.input_path)
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        self.signalp_path = expand_path(self.signalp_path)
        self.signalp3_path = expand_path(self.signalp3_path)
        self.hmmsearch_path = expand_path(self.hmmsearch_path)
        self.blastp_path = expand_path(self.blastp_path)
        self.tmhmm_path = expand_path(self.tmhmm_path)
        self.rxlr_hmm = expand_path(self.rxlr_hmm) if self.rxlr_hmm else DEFAULT_RXLR_HMM
        self.rxlr_wy_hmm = expand_path(self.rxlr_wy_hmm) if self.rxlr_wy_hmm else DEFAULT_RXLR_WY_HMM
        self.crn_hmm = expand_path(self.crn_hmm) if self.crn_hmm else DEFAULT_CRN_HMM
        self.nlp_hmm = expand_path(self.nlp_hmm) if self.nlp_hmm else DEFAULT_NLP_HMM
        self.protease_hmm = expand_path(self.protease_hmm) if self.protease_hmm else DEFAULT_PROTEASE_HMM
        self.scp_hmm = expand_path(self.scp_hmm) if self.scp_hmm else DEFAULT_SCP_HMM
        self.elicitin_hmm = expand_path(self.elicitin_hmm) if self.elicitin_hmm else DEFAULT_ELICITIN_HMM
        self.yxsl_hmm = expand_path(self.yxsl_hmm) if self.yxsl_hmm else DEFAULT_YXSL_HMM
        self.rxlr_blastp_queries = expand_path(self.rxlr_blastp_queries) if self.rxlr_blastp_queries else DEFAULT_RXLR_BLASTP_QUERIES

        # 收集输入文件|Collect input files
        self._collect_input_files()

    def _collect_input_files(self):
        """收集输入FASTA文件|Collect input FASTA files"""
        p = Path(self.input_path)
        fasta_extensions = {'.fa', '.fasta', '.faa', '.pep', '.protein', '.prot'}

        if p.is_file():
            if p.suffix.lower() not in fasta_extensions:
                raise ValueError(f"输入文件应为FASTA格式|Input file should be FASTA format: {p}")
            self._input_files = [str(p)]
        elif p.is_dir():
            self._input_files = sorted([
                str(f) for f in p.iterdir()
                if f.is_file() and f.suffix.lower() in fasta_extensions
            ])
            if not self._input_files:
                raise ValueError(f"目录中未找到FASTA文件|No FASTA files found in directory: {p}")
        else:
            raise ValueError(f"输入路径不存在|Input path not found: {p}")

        # 设置合并FASTA路径|Set combined FASTA path
        if len(self._input_files) == 1:
            self._combined_fasta = self._input_files[0]
        else:
            self._combined_fasta = os.path.join(self.output_dir, '00_pipeline_info', 'combined_input.fasta')
            Path(self._combined_fasta).parent.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        valid_modes = ('rxlr', 'crn', 'nlp', 'protease', 'scp', 'elicitin', 'yxsl')
        if self.mode not in valid_modes:
            errors.append(f"模式必须是{'/'.join(valid_modes)}|Mode must be one of {'/'.join(valid_modes)}: {self.mode}")

        if not self._input_files:
            errors.append("无有效输入文件|No valid input files")

        for f in self._input_files:
            if not os.path.exists(f):
                errors.append(f"输入文件不存在|Input file not found: {f}")

        if not self.skip_signalp and not os.path.exists(self.signalp_path):
            errors.append(f"SignalP程序不存在|SignalP program not found: {self.signalp_path}")

        if self.signalp_version == '3' and not os.path.exists(self.signalp3_path):
            errors.append(f"SignalP 3.0程序不存在|SignalP 3.0 program not found: {self.signalp3_path}")

        if not os.path.exists(self.hmmsearch_path):
            errors.append(f"hmmsearch程序不存在|hmmsearch program not found: {self.hmmsearch_path}")

        if self.rxlr_blastp_queries and not os.path.exists(self.rxlr_blastp_queries):
            errors.append(f"BLASTP查询文件不存在|BLASTP query file not found: {self.rxlr_blastp_queries}")

        if self.rxlr_blastp_queries and not os.path.exists(self.blastp_path):
            errors.append(f"blastp程序不存在|blastp program not found: {self.blastp_path}")

        if self.mode == 'rxlr' and not os.path.exists(self.rxlr_hmm):
            errors.append(f"RxLR HMM文件不存在|RxLR HMM file not found: {self.rxlr_hmm}")

        if self.mode == 'rxlr' and self.use_wy_domain and not os.path.exists(self.rxlr_wy_hmm):
            errors.append(f"RxLR WY HMM文件不存在|RxLR WY HMM file not found: {self.rxlr_wy_hmm}")

        if self.mode == 'rxlr' and not os.path.exists(self.tmhmm_path):
            errors.append(f"tmhmm程序不存在|tmhmm program not found: {self.tmhmm_path}")

        if self.mode == 'crn' and not os.path.exists(self.crn_hmm):
            errors.append(f"CRN HMM文件不存在|CRN HMM file not found: {self.crn_hmm}")

        # 通用效应子类型HMM验证|Generic effector type HMM validation
        type_hmm_map = {
            'nlp': 'nlp_hmm', 'protease': 'protease_hmm',
            'scp': 'scp_hmm', 'elicitin': 'elicitin_hmm', 'yxsl': 'yxsl_hmm',
        }
        if self.mode in type_hmm_map:
            hmm_path = getattr(self, type_hmm_map[self.mode])
            if not os.path.exists(hmm_path):
                errors.append(f"{self.mode.upper()} HMM文件不存在|{self.mode.upper()} HMM file not found: {hmm_path}")

        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
