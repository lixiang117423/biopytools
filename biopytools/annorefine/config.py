"""
annorefine 配置管理|annorefine Configuration Management
注释精修模块的配置类(同源补漏+合并拆分+质控)|Config for annotation refinement
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from ..common.paths import expand_path, get_tool_path, get_samtools_path


@dataclass
class AnnorefineConfig:
    """annorefine 配置类|annorefine Config Class"""

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
    require_complete_orf: bool = True   # partial(覆盖<99)默认丢|drop partial (miniprot coverage based)
    te_overlap_cutoff: float = 50.0     # 真 TE 区重叠阈值 %|TE overlap cutoff
    exclude_te_gap: bool = False        # 质控排除TE区gap(默认不排:真基因可能在TE区)|exclude TE-overlap gaps (default off)

    # ===== 通用生物学质控(普适模块)|General bio-QC =====
    # ① 真实完整 ORF: CDS 长度3倍数 + ATG开头 + 终止密码子结尾(翻译验证, 比 miniprot 覆盖率更严)
    # |Real complete ORF: 3×len + ATG + stop (translation check, stricter than miniprot coverage)
    require_real_orf: bool = True
    # ③ gap 路径坐标零重叠: 与任一BRAKER基因 span 有坐标交集即不算新基因(仅约束纯漏检填补)
    # |Gap path coord-zero-overlap: any coord intersection with a BRAKER gene => not a new gene (gap path only)
    gap_coord_zero_overlap: bool = True
    # ④⑤ 表达证据(用唯一比对 reads 算)|Expression evidence (unique-mapping reads only)
    unique_reads_only: bool = True       # 表达量只用唯一比对 reads(多比对不要, 多为TE/重复假象)|unique reads only
    min_unique_mapq: int = 20            # MAPQ 兜底阈值(samtools 无 -e 时用 -q)|MAPQ fallback
    min_expression_depth: float = 1.0    # 唯一 reads 平均深度下限(>0, 至少有表达)|min unique-read depth (>0)
    min_coverage_breadth: float = 50.0   # CDS 被唯一 reads 覆盖广度 % 下限(防单条reads假象)|min coverage breadth %

    # ===== 路径开关 + 合并拆分判据|Path toggles + merged-gene split =====
    enable_gap_fill: bool = True                   # 纯漏检填补路径(找全新基因)|pure gap-fill path
    enable_split: bool = True                      # 合并拆分路径(拆BRAKER折叠基因)|merged-gene split path
    split_min_hits: int = 2                        # ≥N 个独立命中才判合并
    split_min_copy_coverage: float = 80.0          # 每命中覆盖蛋白≥此%才算完整拷贝

    # ===== 小蛋白回收通道(默认关, 通用)|small-protein lane (default off, general) =====
    # gap_min_cds_len 硬悬崖会把小蛋白整类丢掉; 本通道放宽长度, 用 完整ORF+同源+表达
    # 通用证据找回(无物种/功能假设, 详见 09.二次调试/small_protein_design.md)
    # |Recovers small proteins dropped by gap_min_cds_len, gated by general evidence only
    enable_small_protein: bool = False             # 默认关, 不改现有行为|default off
    small_max_cds_len: int = 450                   # 小蛋白 CDS 上限 bp(150aa)|max CDS (150aa)
    small_min_identity: float = 50.0               # 有表达时放宽 identity|relaxed identity (with expr)
    small_min_coverage: float = 50.0               # 有表达时放宽 coverage|relaxed coverage (with expr)
    small_min_expression_depth: float = 1.0        # 小蛋白表达深度下限|small-protein min depth
    small_min_coverage_breadth: float = 60.0       # 略严于常规 50|stricter breadth
    small_exclude_te: bool = True                  # 强制排 TE|force TE exclusion
    # 强同源直通: identity≥此值的命中(近乎自比对/高度保守)绕过 TE/表达过滤
    # |strong-homology bypass: identity>=this skips TE/expression filters
    # 效应子常在 TE 区且低表达, 强同源(如已知 effector 蛋白自比对)不该被辅助证据拦
    # |effectors often reside in TE regions & low-expression; strong homology
    # (e.g. known effector self-alignment) shouldn't be blocked by auxiliary evidence
    small_strong_homology_identity: float = 95.0   # 默认95(自比对≈100)|default 95 (self-align ~100)
    # 退化模式(无表达数据)同源阈值复用 gap_min_identity/gap_min_coverage(70/80, 更严)
    # |degraded mode (no expression) reuses gap_min_identity/coverage (70/80, stricter)

    # ===== 工具路径(~/... + __post_init__ 展开)|Tool paths =====
    miniprot_bin: str = '~/miniforge3/envs/braker_v.3.0.8/bin/miniprot'
    stringtie_bin: str = '~/.local/bin/stringtie'
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
        self.stringtie_bin = get_tool_path(
            'stringtie', self.stringtie_bin, 'STRINGTIE_PATH')
        if self.samtools_bin:
            self.samtools_bin = expand_path(self.samtools_bin)
        else:
            self.samtools_bin = get_samtools_path()

        # 创建输出目录 + by-step 子目录|Create output + by-step subdirs
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        self.pipeline_info_dir = os.path.join(self.output_dir, '00_pipeline_info')
        self.evidence_dir = os.path.join(self.output_dir, '01_evidence_scan')
        self.gap_dir = os.path.join(self.output_dir, '02_gap_analysis')
        self.gap_filled_dir = os.path.join(self.output_dir, '03_gap_filled')
        self.merged_dir = os.path.join(self.output_dir, '04_merged')
        self.log_dir = os.path.join(self.output_dir, '99_logs')
        for d in [self.pipeline_info_dir, self.evidence_dir, self.gap_dir,
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
