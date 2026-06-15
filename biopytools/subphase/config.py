"""SubPhaser配置管理|SubPhaser Configuration Management"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

from ..common.paths import expand_path

# 亚基因组标签，nsg>26 时自动扩展为 SG1/SG2/...
SUBGENOME_LABELS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


@dataclass
class SubPhaserConfig:
    """SubPhaser亚基因组分离配置|SubPhaser subgenome phasing configuration"""

    # ===== 必需参数|Required parameters =====
    genomes: List[str] = field(default_factory=list)
    nsg: int = 0

    # ===== 输入模式|Input mode =====
    sg_cfgs: Optional[List[str]] = None
    parental_genomes: Optional[List[str]] = None

    # ===== 自动模式标志|Auto mode flag =====
    auto_mode: bool = True

    # ===== 输出配置|Output configuration =====
    output_dir: str = "./subphaser_output"
    prefix: Optional[str] = None

    # ===== 流程参数|Pipeline parameters =====
    threads: int = 24
    overwrite: bool = False
    cleanup: bool = False

    # ===== 染色体过滤|Chromosome filtering =====
    min_chrom_size: int = 1_000_000

    # ===== K-mer参数|K-mer parameters =====
    kmer_size: int = 15
    min_fold: float = 2.0
    min_freq: int = 200

    # ===== 聚类参数|Cluster parameters =====
    max_pval: float = 0.05
    replicates: int = 1000
    test_method: str = "ttest_ind"

    # ===== 步骤控制|Step control =====
    disable_ltr: bool = False
    disable_circos: bool = False
    disable_blocks: bool = False
    just_core: bool = False

    # ===== LTR参数|LTR parameters =====
    ltr_detectors: Optional[List[str]] = None
    mu: float = 13e-9

    # ===== Circos参数|Circos parameters =====
    window_size: int = 1000000
    aligner: str = "minimap2"

    # ===== 高级参数|Advanced parameters =====
    sg_assigned: Optional[str] = None
    target: Optional[str] = None
    labels: Optional[List[str]] = None
    no_label: bool = False
    custom_features: Optional[List[str]] = None
    figfmt: str = "pdf"

    # ===== conda环境|Conda environment =====
    conda_env: str = "SubPhaser"

    # ===== 运行时字段（非参数）|Runtime fields =====
    temp_config: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_dir = expand_path(self.output_dir)
        if self.prefix:
            self.prefix = self.prefix.replace("/", "_")
        if self.sg_cfgs:
            self.sg_cfgs = [expand_path(f) for f in self.sg_cfgs]
        if self.parental_genomes:
            self.parental_genomes = [expand_path(f) for f in self.parental_genomes]
        if self.sg_assigned:
            self.sg_assigned = expand_path(self.sg_assigned)
        if self.target:
            self.target = expand_path(self.target)
        if self.custom_features:
            self.custom_features = [expand_path(f) for f in self.custom_features]

        self.auto_mode = self.sg_cfgs is None

        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not self.genomes:
            errors.append("必须提供基因组文件|Must provide genome files (-i)")
        else:
            for g in self.genomes:
                if not os.path.exists(g):
                    errors.append(f"基因组文件不存在|Genome file not found: {g}")

        if self.nsg < 2:
            errors.append("亚基因组数量必须>=2|Number of subgenomes must be >= 2 (--nsg)")

        if self.sg_cfgs:
            for cfg in self.sg_cfgs:
                if not os.path.exists(cfg):
                    errors.append(f"亚基因组配置文件不存在|Config file not found: {cfg}")

        if self.parental_genomes:
            if len(self.parental_genomes) != 2:
                errors.append(
                    "验证模式需要恰好2个父本基因组|"
                    "Validation mode requires exactly 2 parental genomes"
                )
            else:
                for pg in self.parental_genomes:
                    if not os.path.exists(pg):
                        errors.append(f"父本基因组不存在|Parental genome not found: {pg}")

        if self.sg_assigned and not os.path.exists(self.sg_assigned):
            errors.append(f"亚基因组分配文件不存在|Assignment file not found: {self.sg_assigned}")

        if self.target and not os.path.exists(self.target):
            errors.append(f"目标染色体文件不存在|Target file not found: {self.target}")

        if self.kmer_size < 1:
            errors.append("k-mer大小必须大于0|K-mer size must be > 0")

        if self.min_fold < 1:
            errors.append("最小倍数必须大于1|Min fold must be > 1")

        if self.threads <= 0:
            errors.append("线程数必须为正数|Thread count must be positive")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_subgenome_label(self, index: int) -> str:
        """获取亚基因组标签|Get subgenome label"""
        if index < len(SUBGENOME_LABELS):
            return SUBGENOME_LABELS[index]
        return f"SG{index + 1}"
