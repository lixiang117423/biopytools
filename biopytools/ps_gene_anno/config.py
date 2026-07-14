"""
ps-gene-anno 配置管理|ps-gene-anno Configuration Management
BRAKER 后效应子查漏补缺模块的配置类|Config for post-BRAKER effector gap-filling
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from ..common.paths import expand_path, get_tool_path, get_samtools_path


@dataclass
class PsGeneAnnoConfig:
    """ps-gene-anno 配置类|ps-gene-anno Config Class"""

    # ===== 必需参数|Required =====
    genome: str                 # 未 mask 原始基因组|Unmasked raw genome
    braker_gff3: str             # BRAKER 输出 GFF3|BRAKER output GFF3
    prot_seq: str               # 近缘蛋白|Protein evidence
    output_dir: str             # 输出目录|Output dir

    # ===== 可选证据|Optional evidence =====
    rnaseq_bam: Optional[List[str]] = None    # RNA-seq BAM(辅助深度)|RNA-seq BAMs
    isoseq_bam: Optional[str] = None          # 三代 BAM|Long-read BAM
    repeat_out: Optional[str] = None          # RepeatMasker .out(真TE区排除)|RepeatMasker .out

    # ===== 输出前缀|Output prefix =====
    prefix: Optional[str] = None              # None → genome stem|None => genome stem

    # ===== 质控阈值|QC thresholds =====
    gap_min_identity: float = 70.0      # miniprot identity %|identity cutoff
    gap_min_coverage: float = 80.0      # 命中覆盖蛋白比例 %|coverage cutoff
    gap_min_cds_len: int = 300          # 最小 CDS 长度 bp(过滤短蛋白片段)|min CDS length
    overlap_cutoff: float = 0.0         # 漏检判定:与braker CDS零重叠才算漏检|zero overlap=gap
    require_complete_orf: bool = True   # partial(覆盖<99)默认丢|drop partial
    te_overlap_cutoff: float = 50.0     # 真 TE 区重叠阈值 %|TE overlap cutoff

    # ===== 合并拆分判据|Merged-gene split(保守) =====
    enable_split: bool = True
    split_min_hits: int = 2                        # ≥N 个独立命中才判合并
    split_min_copy_coverage: float = 80.0          # 每命中覆盖蛋白≥此%才算完整拷贝

    # ===== 工具路径(~/... + __post_init__ 展开)|Tool paths =====
    miniprot_bin: str = '~/miniforge3/envs/braker_v.3.0.8/bin/miniprot'
    samtools_bin: str = ''   # 空 → get_samtools_path() 兜底|empty => fallback

    # ===== 流程参数|Pipeline =====
    threads: int = 12

    # ===== 步骤控制|Step control =====
    skip_evidence_scan: bool = False
    skip_gap_analysis: bool = False
    skip_merge: bool = False

    def __post_init__(self):
        """初始化后处理|Post-init: 展开路径 + 创建子目录"""
        # prefix 默认 = genome stem|prefix default = genome stem
        if not self.prefix:
            self.prefix = Path(self.genome).stem

        # 展开用户输入路径(强制绝对路径,避免 CommandRunner 在 output_dir 作 cwd
        # 导致相对路径双重拼接)|Force absolute paths
        def _abs(p):
            return os.path.abspath(os.path.expanduser(os.path.expandvars(p)))
        self.genome = _abs(self.genome)
        self.braker_gff3 = _abs(self.braker_gff3)
        self.prot_seq = _abs(self.prot_seq)
        self.output_dir = _abs(self.output_dir)
        if self.isoseq_bam:
            self.isoseq_bam = _abs(self.isoseq_bam)
        if self.repeat_out:
            self.repeat_out = _abs(self.repeat_out)
        if self.rnaseq_bam:
            self.rnaseq_bam = [_abs(b) for b in self.rnaseq_bam]

        # 展开工具路径|Expand tool paths
        self.miniprot_bin = get_tool_path(
            'miniprot', self.miniprot_bin, 'MINIPROT_PATH')
        if self.samtools_bin:
            self.samtools_bin = expand_path(self.samtools_bin)
        else:
            self.samtools_bin = get_samtools_path()

        # 创建输出目录 + by-step 子目录|Create output + by-step subdirs
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        self.evidence_dir = os.path.join(self.output_dir, '01_evidence_scan')
        self.gap_dir = os.path.join(self.output_dir, '02_gap_analysis')
        self.gap_filled_dir = os.path.join(self.output_dir, '03_gap_filled')
        self.merged_dir = os.path.join(self.output_dir, '04_merged')
        self.log_dir = os.path.join(self.output_dir, '99_logs')
        for d in [self.evidence_dir, self.gap_dir,
                  self.gap_filled_dir, self.merged_dir, self.log_dir]:
            Path(d).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置|Validate configuration"""
        errors = []
        for label, path in [('genome', self.genome),
                            ('braker_gff3', self.braker_gff3),
                            ('prot_seq', self.prot_seq)]:
            if not os.path.exists(path):
                errors.append(f"{label} 不存在|not found: {path}")
        if self.repeat_out and not os.path.exists(self.repeat_out):
            errors.append(f"repeat_out 不存在|not found: {self.repeat_out}")
        if self.threads <= 0:
            errors.append(f"线程数必须为正|threads must be positive: {self.threads}")
        if errors:
            raise ValueError("\n".join(errors))
        return True
