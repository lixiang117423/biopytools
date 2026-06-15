"""
转录组验证注释配置管理模块|Transcriptome Validation Configuration Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Tuple

from ..common.paths import expand_path


@dataclass
class RnaseqValConfig:
    """转录组验证配置类|Transcriptome Validation Configuration Class"""

    # ===== 必需参数|Required parameters =====
    genome_fa: str
    annotation_gtf: str
    output_dir: str = "./rnaseq_val_output"

    # ===== 二代数据|Short-read (2nd gen) data =====
    sr_dir: Optional[str] = None              # 二代 reads 目录（自动检测配对 fastq）
    sr_pattern: Optional[str] = None           # 自定义 fastq 命名模式，如 *_1.clean.fq.gz

    # ===== 三代数据|Long-read (3rd gen) data =====
    lr_dir: Optional[str] = None              # 三代 reads 目录（自动检测单端 fastq）
    lr_platform: str = "pacbio"                # pacbio / ont

    # ===== 通用参数|General parameters =====
    threads: int = 12
    strandness: str = "RF"                    # RF / FR / unstranded
    max_intron: int = 500000
    steps: str = "all"                        # 逗号分隔的步骤名或 all
    sample_timeout: int = 21600               # 单样本超时(秒)，默认6小时

    # ===== StringTie 参数|StringTie parameters =====
    stringtie_min_cov: float = 5.0
    stringtie_min_junction_reads: int = 3
    stringtie_min_isoform_fraction: float = 0.1

    # ===== 校正阈值|Correction thresholds =====
    min_cov: float = 5.0
    min_tpm: float = 0.5
    min_tpm_sr_only: float = 1.0
    min_junction_ont: int = 5

    # ===== 日志选项|Logging options =====
    verbose: bool = False
    quiet: bool = False
    dry_run: bool = False
    force: bool = False
    conda_env: str = "rnaseq_val"              # conda 环境名称|Conda environment name

    # ===== 内部属性|Internal attributes =====
    samples_sr: List[dict] = field(default_factory=list)
    samples_lr: List[dict] = field(default_factory=list)
    output_path: Path = field(default=None, init=False)
    log_dir: Path = field(default=None, init=False)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化输出目录|Normalize output directory
        self.output_dir = os.path.normpath(expand_path(self.output_dir))
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.log_dir = self.output_path / "99_logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)

        # 标准化必需路径|Normalize required paths
        self.genome_fa = os.path.normpath(expand_path(self.genome_fa))
        self.annotation_gtf = os.path.normpath(expand_path(self.annotation_gtf))

        # 标准化可选路径|Normalize optional paths
        if self.sr_dir:
            self.sr_dir = os.path.normpath(expand_path(self.sr_dir))
        if self.lr_dir:
            self.lr_dir = os.path.normpath(expand_path(self.lr_dir))

    def validate(self) -> bool:
        """验证配置参数|Validate configuration parameters

        Returns:
            bool: 验证通过返回 True|True if validation passes

        Raises:
            ValueError: 参数验证失败|Parameter validation failed
        """
        errors = []

        # 检查必需文件|Check required files
        required_files = [
            ('基因组文件|Genome FASTA', self.genome_fa),
            ('注释文件|Annotation GTF', self.annotation_gtf),
        ]
        for desc, path in required_files:
            if not os.path.exists(path):
                errors.append(f"{desc}不存在|does not exist: {path}")

        # 至少需要一种数据源|At least one data source is required
        if not self.sr_dir and not self.lr_dir:
            errors.append(
                "必须提供至少一种数据源(--sr-dir 或 --lr-dir)|"
                "Must provide at least one data source (--sr-dir or --lr-dir)"
            )

        # 检查输入目录|Check input directories
        if self.sr_dir and not os.path.isdir(self.sr_dir):
            errors.append(f"二代数据目录不存在|Short-read directory not found: {self.sr_dir}")
        if self.lr_dir and not os.path.isdir(self.lr_dir):
            errors.append(f"三代数据目录不存在|Long-read directory not found: {self.lr_dir}")

        # 检查 strandness|Check strandness
        valid_strandness = {"RF", "FR", "unstranded", "NONE", "NONE_rf"}
        if self.strandness.upper() not in {s.upper() for s in valid_strandness}:
            errors.append(
                f"无效的链特异性参数|Invalid strandness: {self.strandness}, "
                f"可选|choices: RF, FR, unstranded"
            )

        # 检查 lr_platform|Check lr_platform
        valid_platforms = {"pacbio", "ont"}
        if self.lr_platform not in valid_platforms:
            errors.append(
                f"无效的三代平台参数|Invalid lr_platform: {self.lr_platform}, "
                f"可选|choices: pacbio, ont"
            )

        # 检查线程数|Check threads
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Threads must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_steps_list(self) -> List[str]:
        """解析 steps 参数为步骤列表|Parse steps parameter into a list of step names

        Returns:
            List[str]: 步骤名称列表|List of step names
        """
        all_steps = [
            "align_2nd", "align_3rd",
            "assemble_2nd", "assemble_3rd",
            "compare", "correct", "report",
        ]

        if self.steps.strip().lower() == "all":
            return all_steps

        requested = [s.strip() for s in self.steps.split(",")]

        # 验证步骤名|Validate step names
        invalid = [s for s in requested if s not in all_steps]
        if invalid:
            raise ValueError(
                f"无效的步骤名|Invalid step names: {invalid}\n"
                f"可用步骤|Available steps: {all_steps}"
            )

        return requested

    def get_stringtie_strand_flag(self) -> str:
        """根据 strandness 返回 StringTie 链特异性参数|Return StringTie strand flag based on strandness

        Returns:
            str: StringTie 参数 (--rf, --fr, 或空字符串)|StringTie flag
        """
        s = self.strandness.upper()
        if s == "RF":
            return "--rf"
        elif s == "FR":
            return "--fr"
        return ""

    def get_hisat2_strand_arg(self) -> str:
        """根据 strandness 返回 HISAT2 --rna-strandness 参数|Return HISAT2 strandness argument

        Returns:
            str: HISAT2 参数或空字符串|HISAT2 argument or empty string
        """
        s = self.strandness.upper()
        if s in ("RF", "FR", "NONE", "NONE_RF"):
            return f"--rna-strandness {self.strandness}"
        return ""
