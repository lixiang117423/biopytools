"""
Purge_Dups配置管理模块|Purge_Dups Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
from ..common.paths import expand_path


@dataclass
class PurgeDupsConfig:
    """Purge_Dups配置类|Purge_Dups Configuration Class"""

    # 必需文件|Required files
    input: str  # 基因组组装文件|Genome assembly file
    reads: str  # 可以是PacBio/HiFi reads或Illumina reads列表|Can be PacBio/HiFi reads or Illumina reads list

    # 路径配置|Path configuration
    purge_dups_path: str = '~/miniforge3/envs/purge_dups_v.1.2.6'
    output_dir: str = './purge_dups_output'

    # 处理参数|Processing parameters
    threads: int = 12
    read_type: str = 'hifi'  # 'pacbio', 'hifi', 'illumina'

    # purge_dups参数|purge_dups parameters
    min_fraction: float = 0.8
    min_alignment_score: int = 70
    min_match_score: int = 200
    two_round_chaining: bool = True
    min_match_bases: int = 500
    max_gap_1st: int = 20000
    max_gap_2nd: int = 50000
    min_chain_score: int = 10000
    max_contig_end_ext: int = 15000

    # get_seqs参数|get_seqs parameters
    ends_only: bool = True  # 只去除contig末端的冗余|Only remove duplications at ends
    split_contigs: bool = False
    keep_high_coverage: bool = False
    add_prefix: bool = True
    max_dup_gap: int = 10000
    min_primary_length: int = 10000
    min_length_ratio: float = 0.05

    # calcuts参数|calcuts parameters
    min_depth_fraction: float = 0.1
    manual_cutoffs: Optional[str] = None  # 手动指定阈值文件|Manual cutoffs file

    # 步骤控制|Step control
    step: Optional[int] = None  # 1-5, None表示运行全部步骤|None means run all steps

    # 可选参数|Optional parameters
    split_by_n: bool = False  # split_fa参数: 是否按N分割|split_fa: split by N or not

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开路径中的~和环境变量|Expand ~ and environment variables in paths
        self.purge_dups_path = expand_path(self.purge_dups_path)
        self.input = expand_path(self.input)
        self.output_dir = expand_path(self.output_dir)

        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 创建子目录|Create subdirectories
        self.coverage_dir = self.output_path / "coverage"
        self.split_aln_dir = self.output_path / "split_aln"
        self.purge_dups_dir = self.output_path / "purge_dups"
        self.seqs_dir = self.output_path / "seqs"

        for d in [self.coverage_dir, self.split_aln_dir, self.purge_dups_dir, self.seqs_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.input = os.path.normpath(os.path.abspath(self.input))
        self.purge_dups_path = os.path.normpath(os.path.abspath(self.purge_dups_path))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 设置reads文件|Set reads file(s) - 转换为绝对路径
        if isinstance(self.reads, str):
            self.reads = expand_path(self.reads)
            self.reads = os.path.normpath(os.path.abspath(self.reads))
            self.reads_files = [self.reads]
        else:
            self.reads_files = [
                os.path.normpath(os.path.abspath(expand_path(r))) for r in self.reads
            ]

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.input):
            errors.append(f"基因组文件不存在|Genome file not found: {self.input}")

        for reads_file in self.reads_files:
            if not os.path.exists(reads_file):
                errors.append(f"测序文件不存在|Reads file not found: {reads_file}")

        # 检查purge_dups路径|Check purge_dups path
        if not os.path.exists(self.purge_dups_path):
            errors.append(f"Purge_Dups路径不存在|Purge_Dups path does not exist: {self.purge_dups_path}")

        # 检查read_type|Check read_type
        valid_read_types = ['pacbio', 'hifi', 'illumina']
        if self.read_type not in valid_read_types:
            errors.append(f"无效的reads类型|Invalid reads type: {self.read_type} (应为|should be one of {valid_read_types})")

        # 检查步骤参数|Check step parameter
        if self.step is not None and self.step not in [1, 2, 3, 4, 5]:
            errors.append(f"无效的步骤编号|Invalid step number: {self.step} (应为1-5|should be 1-5)")

        # 检查参数范围|Check parameter ranges
        if not 0 <= self.min_fraction <= 1:
            errors.append(f"min_fraction必须在0-1之间|min_fraction must be between 0 and 1")

        if not 0 <= self.min_length_ratio <= 1:
            errors.append(f"min_length_ratio必须在0-1之间|min_length_ratio must be between 0 and 1")

        if errors:
            raise ValueError("\n".join(errors))

        return True
