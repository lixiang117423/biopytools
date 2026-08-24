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
    clean_fastq_dir: Optional[str] = None   # 已 clean 的 fastq(给则跳过 QC)|clean fastq (skip QC)
    genome: str = ""
    output_dir: str = "mixrace_out"
    repeat_bed: Optional[str] = None     # L4 额外排除区域BED(与自动热点并集)|extra exclude regions
    # 寄主剔除|host depletion
    host_genome: Optional[str] = None    # 给则 clean 后比对寄主并整对剔除|deplete host reads if set
    min_mapq: int = 20                   # 「比对上」判定阈值+mapped reads提取过滤|MAPQ threshold
    # 三分支判读阈值(v0.3)|three-branch verdict thresholds
    pure_het_threshold: float = 0.001    # 总杂合率<0.1% 判纯;同作 L4 热点候选筛选口径(双语义)
                                         # |pure verdict threshold; also the L4 hotspot-candidate filter
    partner_alt_min: float = 0.8         # 伴侣携带ALT率阈值|partner ALT-carrier min
    partner_hom_min: float = 0.5         # 伴侣纯合1/1占比阈值|partner homozygous min
    min_sites: int = 1000                # 最低有GT位点数|min called sites
    # 热点|hotspots
    window_size: int = 100000
    hotspot_fold: float = 2.0
    hotspot_min_median: float = 0.10
    # 执行|Execution
    threads: int = 12
    sample_parallel: int = 1        # 样本级并行数(每worker线程=threads/该值)|per-sample workers
    kmer_size: int = 21
    read_length: int = 150
    step: Optional[int] = None        # None=全跑|None=all steps, 1..5=单步|single step
    enable_checkpoint: bool = True
    dry_run: bool = False
    # 工具路径(默认~, __post_init__ 展开)|Tool paths (default ~, expanded in __post_init__)
    bwa_mem2_path: str = "~/miniforge3/envs/cphasing/bin/bwa-mem2"
    samtools_path: str = "~/miniforge3/envs/align/bin/samtools"
    bcftools_path: str = "~/miniforge3/envs/align/bin/bcftools"

    def __post_init__(self):
        # 展开工具路径(支持环境变量覆盖)|expand tool paths (env var override supported)
        self.bwa_mem2_path = get_tool_path("bwa-mem2", self.bwa_mem2_path, "BWA_MEM2_PATH")
        self.samtools_path = get_tool_path("samtools", self.samtools_path, "SAMTOOLS_PATH")
        self.bcftools_path = get_tool_path("bcftools", self.bcftools_path, "BCFTOOLS_PATH")
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
        if self.host_genome:
            self.host_genome = expand_path(self.host_genome)

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
        if self.host_genome and not os.path.isfile(self.host_genome):
            errors.append(f"寄主基因组不存在|Host genome not found: {self.host_genome}")
        if self.min_mapq < 0:
            errors.append("min_mapq不能为负|min_mapq must be >= 0")
        if self.pure_het_threshold <= 0 or self.pure_het_threshold >= 1:
            errors.append("pure_het_threshold须在(0,1)|pure_het_threshold must be in (0,1)")
        if self.min_sites <= 0:
            errors.append("min_sites必须为正|min_sites must be positive")
        if self.sample_parallel < 1:
            errors.append("sample_parallel须>=1|sample_parallel must be >= 1")
        if self.window_size <= 0:
            errors.append("window_size必须为正|window_size must be positive")
        if not 0 <= self.partner_alt_min <= 1:
            errors.append("partner_alt_min须在[0,1]|partner_alt_min must be in [0,1]")
        if not 0 <= self.partner_hom_min <= 1:
            errors.append("partner_hom_min须在[0,1]|partner_hom_min must be in [0,1]")
        if self.hotspot_fold < 1:
            errors.append("hotspot_fold必须>=1|hotspot_fold must be >= 1")
        if not 0 <= self.hotspot_min_median <= 1:
            errors.append("hotspot_min_median须在[0,1]|hotspot_min_median must be in [0,1]")
        if self.step is not None and not (1 <= self.step <= 5):
            errors.append("step须为1-5|step must be 1-5")
        if errors:
            raise ValueError("\n".join(errors))
