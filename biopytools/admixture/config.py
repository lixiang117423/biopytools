"""
ADMIXTURE分析配置管理模块|ADMIXTURE Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
from ..common.paths import expand_path


@dataclass
class AdmixtureConfig:
    """ADMIXTURE分析配置类|ADMIXTURE Analysis Configuration Class"""

    # 必需文件|Required files
    vcf_file: str
    output_dir: str = "admixture_results"

    # 分析方法|Analysis method
    method: str = "admixture"  # "admixture" or "adamixture"

    # 分析参数|Analysis parameters
    min_k: int = 2
    max_k: int = 10
    cv_folds: int = 5
    threads: int = 64

    # ADAMIXTURE 参数|ADAMIXTURE parameters
    adamixture_path: str = "~/miniforge3/envs/adamixture_v.1.0.2/bin/adamixture"
    adamixture_lr: float = 0.005
    adamixture_beta1: float = 0.80
    adamixture_beta2: float = 0.88
    adamixture_max_iter: int = 1500
    adamixture_seed: int = 42

    # 质控参数|Quality control parameters
    maf: float = 0.05
    missing_rate: float = 0.1
    hwe_pvalue: float = 1e-6

    # LD剪枝参数|LD pruning parameters
    ld_prune: bool = True
    ld_window: str = "3000kb"   # 窗口，支持kb或SNP数|window in kb or SNP count
    ld_step: int = 1            # 步长|step size
    ld_r2: float = 0.2          # r2阈值|r2 threshold

    # 处理选项|Processing options
    skip_preprocessing: bool = False
    keep_intermediate: bool = False

    # 日志配置|Logging configuration
    log_level: str = "INFO"
    quiet: bool = False
    verbose: int = 0

    # 执行控制|Execution control
    force: bool = False
    dry_run: bool = False

    # 内部属性|Internal attributes
    base_name: str = "admixture_ready"

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.adamixture_path = expand_path(self.adamixture_path)
        self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在|VCF file does not exist: {self.vcf_file}")

        # 检查方法参数|Check method parameter
        if self.method not in ["admixture", "adamixture"]:
            errors.append(f"方法必须是 'admixture' 或 'adamixture'|Method must be 'admixture' or 'adamixture': {self.method}")

        # 检查K值范围|Check K range
        if self.min_k < 1 or self.max_k < self.min_k:
            errors.append(f"无效的K值范围|Invalid K range: {self.min_k} to {self.max_k}")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive: {self.threads}")

        # 检查质控参数|Check QC parameters
        if not 0 <= self.maf <= 0.5:
            errors.append(f"MAF值应在0-0.5之间|MAF should be between 0-0.5: {self.maf}")

        if not 0 <= self.missing_rate <= 1:
            errors.append(f"缺失率应在0-1之间|Missing rate should be between 0-1: {self.missing_rate}")

        # 检查LD剪枝参数|Check LD pruning parameters
        if not 0 < self.ld_r2 < 1:
            errors.append(f"LD剪枝r2阈值应在0-1之间|LD pruning r2 should be between 0 and 1: {self.ld_r2}")

        if errors:
            raise ValueError("配置错误|Configuration error:\n" + "\n".join(errors))

        return True
