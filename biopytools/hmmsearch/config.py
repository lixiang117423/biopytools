"""
HMMsearch配置管理模块|HMMsearch Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional
from ..common.paths import expand_path

FASTA_EXTENSIONS = ('.fa', '.faa', '.fasta', '.fa.gz', '.faa.gz', '.fasta.gz')


@dataclass
class HMMsearchConfig:
    """HMMsearch配置类|HMMsearch Configuration Class"""

    # 必需参数|Required parameters (模式1: 处理已有domtblout)
    domtblout_file: str = None
    protein_fasta: str = None

    # 必需参数|Required parameters (模式2: 运行hmmsearch)
    hmm_file: str = None

    # 输出配置|Output configuration
    output_dir: str = './hmmsearch_output'
    output_prefix: str = 'hmmsearch_results'

    # HMMsearch软件配置|HMMsearch software configuration
    hmmsearch_path: str = '~/miniforge3/envs/resistify_v.1.3.0/bin/hmmsearch'
    threads: int = 12

    # HMMsearch运行参数|HMMsearch run parameters
    evalue_cutoff: float = None
    score_cutoff: float = None
    use_cut_tc: bool = False
    use_cut_ga: bool = False
    use_cut_nc: bool = False

    # 过滤参数|Filtering parameters (用于结果过滤，区别于hmmsearch运行参数)
    evalue_threshold: float = None
    score_threshold: float = None

    # 序列提取选项|Sequence extraction options
    extract_protein_sequences: bool = True
    extract_domain_sequences: bool = True

    # 是否输出CSV和Excel|Whether to output CSV and Excel
    output_csv: bool = True
    output_excel: bool = True

    # 展开后的蛋白序列文件列表（支持单个文件或目录）|Expanded protein fasta file list (supports single file or directory)
    protein_fastas: List[str] = field(default_factory=list)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        if self.domtblout_file:
            self.domtblout_file = os.path.normpath(os.path.abspath(self.domtblout_file))
        if self.protein_fasta:
            self.protein_fasta = os.path.normpath(os.path.abspath(self.protein_fasta))
            if os.path.isdir(self.protein_fasta):
                self.protein_fastas = self._collect_fasta_files(self.protein_fasta)
            else:
                self.protein_fastas = [self.protein_fasta]
        if self.hmm_file:
            self.hmm_file = os.path.normpath(os.path.abspath(self.hmm_file))
        if self.hmmsearch_path:
            self.hmmsearch_path = os.path.normpath(os.path.abspath(expand_path(self.hmmsearch_path)))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

    def _collect_fasta_files(self, directory: str) -> List[str]:
        """收集目录中的FASTA文件|Collect FASTA files from directory"""
        fasta_files = []
        for f in sorted(os.listdir(directory)):
            if any(f.lower().endswith(ext) for ext in FASTA_EXTENSIONS):
                fasta_files.append(os.path.normpath(os.path.join(directory, f)))
        return fasta_files

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查模式|Check mode
        mode1 = self.domtblout_file is not None
        mode2 = self.hmm_file is not None and self.protein_fasta is not None

        if not mode1 and not mode2:
            errors.append("必须指定domtblout文件，或者同时指定hmm文件和蛋白序列文件|Must specify domtblout file, or both hmm file and protein fasta file")
        elif mode1 and mode2:
            errors.append("不能同时指定domtblout文件和hmm文件|Cannot specify both domtblout file and hmm file")

        # 模式1验证|Mode 1 validation (process existing domtblout)
        if mode1:
            if not os.path.exists(self.domtblout_file):
                errors.append(f"domtblout文件不存在|domtblout file not found: {self.domtblout_file}")
            if self.protein_fasta and not os.path.exists(self.protein_fasta):
                errors.append(f"蛋白序列路径不存在|Protein FASTA path not found: {self.protein_fasta}")

        # 模式2验证|Mode 2 validation (run hmmsearch)
        if mode2:
            if not os.path.exists(self.hmm_file):
                errors.append(f"HMM文件不存在|HMM file not found: {self.hmm_file}")
            if not os.path.exists(self.protein_fasta):
                errors.append(f"蛋白序列路径不存在|Protein FASTA path not found: {self.protein_fasta}")
            if os.path.isdir(self.protein_fasta) and not self.protein_fastas:
                errors.append(f"目录中未找到FASTA文件|No FASTA files found in directory: {self.protein_fasta}")
            if not os.path.exists(self.hmmsearch_path):
                errors.append(f"hmmsearch程序不存在|hmmsearch program not found: {self.hmmsearch_path}")

        # 检查阈值参数|Check threshold parameters
        if self.evalue_threshold is not None and self.evalue_threshold <= 0:
            errors.append(f"E-value阈值必须为正数|E-value threshold must be positive: {self.evalue_threshold}")

        if self.score_threshold is not None and self.score_threshold < 0:
            errors.append(f"分数阈值不能为负数|Score threshold cannot be negative: {self.score_threshold}")

        # 检查cut选项互斥|Check cut options mutual exclusion
        cut_options = sum([self.use_cut_tc, self.use_cut_ga, self.use_cut_nc])
        if cut_options > 1:
            errors.append("cut_tc、cut_ga、cut_nc不能同时使用|cut_tc, cut_ga, cut_nc cannot be used simultaneously")

        if errors:
            raise ValueError("\n".join(errors))

        return True
