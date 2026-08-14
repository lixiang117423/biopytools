"""mixrace 配置|mixrace configuration

所有 ~ 路径在 __post_init__ 中经 get_tool_path 展开。
|All ~ paths are expanded via get_tool_path in __post_init__.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional
import os

from ..common.paths import expand_path, get_tool_path


@dataclass
class MixraceConfig:
    """mixrace 全流程配置|mixrace pipeline configuration."""

    # 输入|Input
    fastq_dir: str = ""
    clean_fastq_dir: Optional[str] = None   # 已 clean 的 fastq(给则跳过 QC,freebayes calling 与 smudgescope 都基于它)
    genome: str = ""
    output_dir: str = "mixrace_out"
    repeat_bed: Optional[str] = None
    # 执行|Execution
    threads: int = 12
    kmer_size: int = 21
    read_length: int = 150
    step: Optional[int] = None        # None=全跑|None=all steps, 1..7=单步|single step
    enable_checkpoint: bool = True
    dry_run: bool = False
    # 工具路径(默认~, __post_init__ 展开)|Tool paths (default ~, expanded in __post_init__)
    bwa_mem2_path: str = "~/miniforge3/envs/cphasing/bin/bwa-mem2"
    samtools_path: str = "~/miniforge3/envs/sv_calling/bin/samtools"
    bcftools_path: str = "~/miniforge3/envs/sv_calling/bin/bcftools"
    bedtools_path: str = "~/miniforge3/envs/sv_calling/bin/bedtools"
    rscript_path: str = "~/miniforge3/envs/WGCNA_v.1.73/bin/Rscript"
    freebayes_path: str = "~/miniforge3/envs/freebayes/bin/freebayes"
    # freebayes 单倍体 calling 参数(导师: -p 1 --min-alternate-fraction 0.02 --min-coverage 30)
    freebayes_min_alternate_fraction: float = 0.02
    freebayes_min_coverage: int = 30
    # 过滤|filtering
    min_qual: int = 30
    min_dp: int = 15
    min_alt_reads: int = 3
    # 杂合率判读阈值(单倍体,导师)|het-rate verdict thresholds (haploid, advisor)
    het_pure: float = 0.0001         # <0.01% 纯|pure
    het_suspicious: float = 0.001    # 0.1% 可疑|suspicious
    het_impure: float = 0.01         # 1% 不纯|impure
    min_depth: float = 50.0          # 导师:严格过滤需 ≥50x|advisor: >=50x for strict filtering
    # 校准(已知纯样品)|calibration (known-pure samples)
    pure_samples: Optional[List[str]] = None

    def __post_init__(self):
        # 展开工具路径(支持环境变量覆盖)|expand tool paths (env var override supported)
        self.bwa_mem2_path = get_tool_path("bwa-mem2", self.bwa_mem2_path, "BWA_MEM2_PATH")
        self.freebayes_path = get_tool_path("freebayes", self.freebayes_path, "FREEBAYES_PATH")
        self.samtools_path = get_tool_path("samtools", self.samtools_path, "SAMTOOLS_PATH")
        self.bcftools_path = get_tool_path("bcftools", self.bcftools_path, "BCFTOOLS_PATH")
        self.bedtools_path = get_tool_path("bedtools", self.bedtools_path, "BEDTOOLS_PATH")
        self.rscript_path = get_tool_path("Rscript", self.rscript_path, "RSCRIPT_PATH")
        # 用户输入路径统一用 expand_path(同时展开 ~ 和 $VAR;裸 expanduser 漏 $VAR)
        # |User paths via expand_path (expands ~ AND $VAR; bare expanduser misses $VAR)
        if self.fastq_dir:
            self.fastq_dir = expand_path(self.fastq_dir)
        if self.clean_fastq_dir:
            self.clean_fastq_dir = expand_path(self.clean_fastq_dir)
        if self.genome:
            self.genome = expand_path(self.genome)
        if self.output_dir:
            self.output_dir = expand_path(self.output_dir)
        if self.repeat_bed:
            self.repeat_bed = expand_path(self.repeat_bed)

    def validate(self):
        """校验配置,错误一次性抛出|validate config, raise all errors at once."""
        errors = []
        # 输入二选一:raw fastq_dir 或 clean_fastq_dir|one of raw fastq_dir / clean_fastq_dir
        if self.clean_fastq_dir:
            if not os.path.isdir(self.clean_fastq_dir):
                errors.append(f"clean fastq目录不存在|clean fastq dir not found: {self.clean_fastq_dir}")
        elif not self.fastq_dir or not os.path.isdir(self.fastq_dir):
            errors.append(f"输入fastq目录不存在|Input fastq dir not found: {self.fastq_dir}")
        if not self.genome or not os.path.isfile(self.genome):
            errors.append(f"参考基因组不存在|Reference genome not found: {self.genome}")
        if self.threads <= 0:
            errors.append("线程数必须为正|Threads must be positive")
        if self.repeat_bed and not os.path.isfile(self.repeat_bed):
            errors.append(f"repeat bed不存在|Repeat bed not found: {self.repeat_bed}")
        if self.step is not None and not (1 <= self.step <= 7):
            errors.append("step须为1-7|step must be 1-7")
        if errors:
            raise ValueError("\n".join(errors))
