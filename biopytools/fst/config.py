"""
Fst计算配置管理模块|Fst Calculation Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path, get_tool_path


@dataclass
class FstConfig:
    """Fst计算配置类|Fst Calculation Configuration Class"""

    # 必需参数|Required parameters
    vcf_file: str
    pop_file: str
    output_dir: str

    # 软件路径|Software paths
    plink_path: str = None

    # 质控参数|Quality control parameters
    enable_qc: bool = False  # 是否进行质控过滤|Whether to perform quality control filtering
    maf: float = 0.05
    geno: float = 0.1
    mind: float = 0.1
    hwe: float = 1e-6

    # 输出控制|Output control
    keep_intermediate: bool = True

    # Bootstrap参数|Bootstrap parameters
    enable_bootstrap: bool = False
    bootstrap_iterations: int = 100
    min_samples: int = 10  # 排除样本数少于此值的群体|Exclude populations with fewer samples
    exclude_pops: Optional[str] = None  # 手动指定要排除的群体（逗号分隔）|Manually specify populations to exclude (comma-separated)

    # LD pruning参数|LD pruning parameters
    enable_ld_prune: bool = True
    ld_window_size: int = 50
    ld_step_size: int = 10
    ld_r2_threshold: float = 0.2

    # SNP抽稀参数|SNP thinning parameters
    thin_threshold: Optional[float] = None  # SNP抽稀比例（如0.5表示保留50%）|SNP thinning ratio (e.g., 0.5 means keep 50%)

    # 并行参数|Parallel parameters
    threads: int = 12

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 设置plink默认路径|Set default plink path
        if self.plink_path is None:
            self.plink_path = get_tool_path(
                'plink',
                '~/miniforge3/envs/Population_genetics/bin/plink',
                'PLINK_PATH'
            )

        # 展开路径|Expand paths
        self.vcf_file = expand_path(self.vcf_file)
        self.pop_file = expand_path(self.pop_file)
        self.output_dir = expand_path(self.output_dir)
        self.plink_path = expand_path(self.plink_path)

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir).resolve()
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 设置中间文件目录|Set intermediate files directory
        self.intermediate_dir = self.output_path / '00_intermediate'
        if self.keep_intermediate:
            self.intermediate_dir.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在|VCF file does not exist: {self.vcf_file}")

        if not os.path.exists(self.pop_file):
            errors.append(f"群体文件不存在|Population file does not exist: {self.pop_file}")

        # 检查plink路径|Check plink path
        if not os.path.exists(self.plink_path):
            errors.append(f"PLINK路径不存在|PLINK path does not exist: {self.plink_path}")

        # 检查质控参数范围|Check QC parameter ranges
        if not 0 <= self.maf <= 1:
            errors.append(f"MAF阈值必须在0-1之间|MAF threshold must be between 0-1: {self.maf}")

        if not 0 <= self.geno <= 1:
            errors.append(f"GENO阈值必须在0-1之间|GENO threshold must be between 0-1: {self.geno}")

        if not 0 <= self.mind <= 1:
            errors.append(f"MIND阈值必须在0-1之间|MIND threshold must be between 0-1: {self.mind}")

        if not 0 <= self.hwe <= 1:
            errors.append(f"HWE阈值必须在0-1之间|HWE threshold must be between 0-1: {self.hwe}")

        # 检查bootstrap参数|Check bootstrap parameters
        if self.bootstrap_iterations < 1:
            errors.append(f"Bootstrap迭代次数必须>=1|Bootstrap iterations must be >= 1: {self.bootstrap_iterations}")

        if self.min_samples < 2:
            errors.append(f"最小样本数必须>=2|Minimum samples must be >= 2: {self.min_samples}")

        # 检查LD pruning参数|Check LD pruning parameters
        if not 0 < self.ld_r2_threshold <= 1:
            errors.append(f"LD R2阈值必须在0-1之间|LD R2 threshold must be between 0-1: {self.ld_r2_threshold}")

        # 检查thin阈值|Check thin threshold
        if self.thin_threshold is not None and not 0 < self.thin_threshold <= 1:
            errors.append(f"SNP抽稀阈值必须在0-1之间|SNP thinning threshold must be between 0-1: {self.thin_threshold}")

        # 检查线程数|Check threads
        if self.threads < 1:
            errors.append(f"线程数必须>=1|Threads must be >= 1: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
