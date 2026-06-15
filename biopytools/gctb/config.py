"""
GCTB配置管理模块|GCTB Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path, get_tool_path

@dataclass
class GCTBConfig:
    """GCTB配置类|GCTB Configuration Class"""

    # 必需参数|Required parameters
    vcf_file: str
    pheno_file: str

    # 输出配置|Output configuration
    output_dir: str = "./gctb_output"

    # 软件路径配置|Software path configuration
    gctb_path: str = field(
        default_factory=lambda: get_tool_path(
            'gctb',
            '~/miniforge3/envs/gctb/bin/gctb',
            'GCTB_PATH'
        )
    )
    plink_path: str = field(
        default_factory=lambda: get_tool_path(
            'plink',
            '~/miniforge3/envs/Population_genetics/bin/plink',
            'PLINK_PATH'
        )
    )

    # 质量控制参数|Quality control parameters
    maf_threshold: float = 0.01  # MAF阈值|MAF threshold
    miss_threshold: float = 0.1  # 缺失率阈值|Missing rate threshold
    hwe_p: float = 1e-6  # HWE p值阈值|HWE p-value threshold

    # 分析参数|Analysis parameters
    bayes_type: str = "S"  # 贝叶斯模型类型: S, R, C|Bayesian model type
    analysis_mode: str = "individual"  # 分析模式: individual, summary|Analysis mode
    ld_matrix_type: str = "sparse"  # LD矩阵类型: sparse, block, eigen|LD matrix type

    # GCTB 高级参数|GCTB advanced parameters
    threads: int = 12  # 线程数|Number of threads
    seed: Optional[int] = None  # 随机种子|Random seed
    pi: Optional[float] = None  # polygenicity参数|Polygenicity parameter
    sigma_g: Optional[float] = None  # 遗传方差|Genetic variance
    rho: Optional[str] = None  # SNP效应与MAF关系参数|Parameter for effect-MAF relationship

    # 步骤控制|Step control
    step: Optional[str] = None  # 运行指定步骤: convert, qc, ld, analysis|Run specific step only
    batch: bool = True  # 批量处理多个表型（默认开启）|Batch process multiple phenotypes (enabled by default)
    keep_intermediate: bool = False  # 保留中间文件|Keep intermediate files

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有路径|Expand all paths
        self.vcf_file = expand_path(self.vcf_file)
        self.pheno_file = expand_path(self.pheno_file)
        self.output_dir = expand_path(self.output_dir)
        self.gctb_path = expand_path(self.gctb_path)
        self.plink_path = expand_path(self.plink_path)

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 创建子目录|Create subdirectories
        self.convert_dir = self.output_path / "01_vcf_to_plink"
        self.qc_dir = self.output_path / "02_quality_control"
        self.ld_dir = self.output_path / "03_ld_matrix"
        self.analysis_dir = self.output_path / "04_gctb_analysis"
        self.logs_dir = self.output_path / "99_logs"

        for dir_path in [self.convert_dir, self.qc_dir, self.ld_dir,
                         self.analysis_dir, self.logs_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在|VCF file not found: {self.vcf_file}")
        if not os.path.exists(self.pheno_file):
            errors.append(f"表型文件不存在|Phenotype file not found: {self.pheno_file}")

        # 检查软件路径|Check software paths
        if not os.path.exists(self.gctb_path):
            errors.append(f"GCTB不存在|GCTB not found: {self.gctb_path}")
        if not os.path.exists(self.plink_path):
            errors.append(f"PLINK不存在|PLINK not found: {self.plink_path}")

        # 验证参数范围|Validate parameter ranges
        if not 0 < self.maf_threshold <= 0.5:
            errors.append(f"MAF阈值必须在0-0.5之间|MAF threshold must be between 0 and 0.5: {self.maf_threshold}")
        if not 0 < self.miss_threshold <= 1:
            errors.append(f"缺失率阈值必须在0-1之间|Missing rate threshold must be between 0 and 1: {self.miss_threshold}")

        # 验证贝叶斯模型类型|Validate Bayesian model type
        if self.bayes_type.upper() not in ['S', 'R', 'C']:
            errors.append(f"贝叶斯模型类型必须是S、R或C|Bayesian model type must be S, R, or C: {self.bayes_type}")

        # 验证分析模式|Validate analysis mode
        if self.analysis_mode not in ['individual', 'summary']:
            errors.append(f"分析模式必须是individual或summary|Analysis mode must be individual or summary: {self.analysis_mode}")

        # 验证LD矩阵类型|Validate LD matrix type
        if self.ld_matrix_type not in ['sparse', 'block', 'eigen']:
            errors.append(f"LD矩阵类型必须是sparse、block或eigen|LD matrix type must be sparse, block, or eigen: {self.ld_matrix_type}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
