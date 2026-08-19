"""easyhap 配置|easyhap configuration"""
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path, get_tool_path

REGION_RE = re.compile(r"^([^:]+):(\d+)-(\d+)$")

MODES = ("inbred", "hybrid")
HETERO_POLICIES = ("slash", "iupac", "missing")
FISHER_ADJUSTS = ("none", "bh")
VCF_BACKENDS = ("auto", "cyvcf2", "pysam", "plain")
PLOT_HAP_LEVELS = ("hap", "cluster")
PLOT_FORMATS = ("pdf", "svg", "png")

# 存在性检查清单: (字段名, 中文名|英文名)|existence check list
FILE_LABELS = (
    ("input_vcf", "输入VCF|Input VCF"),
    ("group_file", "分组表|Group file"),
    ("region_file", "区域文件|Region file"),
    ("traits", "性状表|Trait file"),
    ("gff", "GFF注释|GFF annotation"),
)


@dataclass
class EasyHapConfig:
    """EasyHap 区域单倍型分析配置|EasyHap regional haplotype config"""

    # 必需|Required
    input_vcf: str
    group_file: str
    region: Optional[str] = None        # 单区域 chr:start-end|single region
    region_file: Optional[str] = None   # 批量区域文件|batch region file

    # 输出|Output
    output_dir: str = "./easyhap_output"

    # 分析|Analysis
    mode: str = "inbred"
    hetero_policy: str = "slash"
    cluster_threshold: float = 0.15
    vcf_backend: str = "auto"
    no_processed: bool = False

    # 过滤|Filtering
    fisher_groups: Optional[str] = None
    fisher_alpha: Optional[float] = None
    fisher_adjust: str = "none"

    # 绘图|Plotting
    plot: bool = False
    gff: Optional[str] = None
    traits: Optional[str] = None
    trait_cols: Optional[str] = None
    plot_format: str = "pdf"
    plot_hap_level: str = "hap"
    plot_min_count: int = 1

    # 运行|Runtime
    force: bool = False
    log_level: str = "INFO"
    log_file: Optional[str] = None

    # 工具路径|Tool path
    easyhap_path: str = field(
        default_factory=lambda: get_tool_path(
            "easyhap", "~/miniforge3/envs/pop/bin/easyhap", "EASYHAP_PATH"))

    def __post_init__(self):
        """展开路径、建目录|Expand paths, create dirs"""
        self.input_vcf = expand_path(self.input_vcf)
        self.group_file = expand_path(self.group_file)
        self.output_dir = expand_path(self.output_dir)
        for attr in ("region_file", "traits", "gff"):
            value = getattr(self, attr)
            if value:
                setattr(self, attr, expand_path(value))
        if self.log_file:
            self.log_file = expand_path(self.log_file)

        self.output_path = Path(self.output_dir)
        self.haplotypes_dir = self.output_path / "01_haplotypes"
        self.info_dir = self.output_path / "00_pipeline_info"
        self.logs_dir = self.output_path / "99_logs"
        for d in (self.haplotypes_dir, self.info_dir, self.logs_dir):
            d.mkdir(parents=True, exist_ok=True)
        if not self.log_file:
            self.log_file = str(self.logs_dir / "easyhap.log")

    def validate(self) -> bool:
        """校验配置(错误一次性收集)|Validate; raise ValueError with all errors"""
        errors = []
        for attr, label in FILE_LABELS:
            path = getattr(self, attr)
            if path and not Path(path).exists():
                errors.append(f"{label}不存在|not found: {path}")
        if not self.region and not self.region_file:
            errors.append("必须提供 --region 或 --region-file|"
                          "Must provide --region or --region-file")
        if self.region and self.region_file:
            errors.append("--region 与 --region-file 互斥|"
                          "--region and --region-file are mutually exclusive")
        if self.region:
            m = REGION_RE.match(self.region.strip())
            if not m:
                errors.append(f"区域格式错误(需 chr:start-end)|"
                              f"Bad region format (need chr:start-end): {self.region!r}")
            elif int(m.group(2)) > int(m.group(3)):
                errors.append(f"区域起始大于终止|Region start > end: {self.region!r}")
        if self.mode not in MODES:
            errors.append(f"无效模式|Invalid mode: {self.mode} (inbred/hybrid)")
        if self.hetero_policy not in HETERO_POLICIES:
            errors.append(f"无效杂合策略|Invalid hetero policy: {self.hetero_policy}")
        if self.vcf_backend not in VCF_BACKENDS:
            errors.append(f"无效后端|Invalid backend: {self.vcf_backend}")
        if self.fisher_adjust not in FISHER_ADJUSTS:
            errors.append(f"无效校正方法|Invalid fisher adjust: {self.fisher_adjust}")
        if self.plot_hap_level not in PLOT_HAP_LEVELS:
            errors.append(f"无效绘图层级|Invalid plot level: {self.plot_hap_level}")
        if self.fisher_groups:
            parts = [x.strip() for x in self.fisher_groups.split(",") if x.strip()]
            if len(parts) != 2:
                errors.append("--fisher-groups 需恰两个逗号分隔的分组名|"
                              "--fisher-groups needs exactly two comma-separated names")
        if self.fisher_alpha is not None and not (0 < self.fisher_alpha <= 1):
            errors.append(f"--fisher-alpha 须在 (0,1]|"
                          f"--fisher-alpha must be in (0,1]: {self.fisher_alpha}")
        if not (0 <= self.cluster_threshold <= 1):
            errors.append(f"--cluster-threshold 须在 [0,1]|"
                          f"--cluster-threshold must be in [0,1]: {self.cluster_threshold}")
        if self.plot_min_count < 1:
            errors.append(f"--plot-min-count 须 >= 1|"
                          f"--plot-min-count must be >= 1: {self.plot_min_count}")
        if self.trait_cols and not self.traits:
            errors.append("--trait-cols 需要 --traits|--trait-cols requires --traits")
        for fmt in [f.strip() for f in self.plot_format.split(",") if f.strip()]:
            if fmt not in PLOT_FORMATS:
                errors.append(f"无效图格式|Invalid plot format: {fmt} (pdf/svg/png)")
        if errors:
            raise ValueError("\n".join(errors))
        return True
