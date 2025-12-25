"""
📋 GWAS分析配置管理模块 | GWAS Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple


@dataclass
class GWASConfig:
    """GWAS分析配置类 | GWAS Analysis Configuration Class"""

    # 输入文件 | Input files
    vcf_file: Optional[str] = None
    phenotype_file: Optional[str] = None

    # 输出文件 | Output files
    output_prefix: str = "gwas_analysis"
    manhattan_plot: Optional[str] = None
    qq_plot: Optional[str] = None
    significant_snps: Optional[str] = None
    excel_report: Optional[str] = None
    html_report: Optional[str] = None

    # 工具路径 | Tool paths
    plink_path: str = "plink"
    bcftools_path: str = "bcftools"
    admixture_path: str = "admixture"
    emmax_path: str = "emmax"

    # VCF过滤参数 | VCF filtering parameters
    maf_threshold: float = 0.05  # 最小等位基因频率 | Minor allele frequency threshold
    missing_threshold: float = 0.2  # 缺失率阈值 | Missing rate threshold
    depth_min: int = 3  # 最小测序深度 | Minimum sequencing depth
    depth_max: int = 50  # 最大测序深度 | Maximum sequencing depth
    qual_min: float = 20.0  # 最小质量值 | Minimum quality value

    # 群体结构分析参数 | Population structure analysis parameters
    admixture_k_range: Tuple[int, int] = (
        1,
        20,
    )  # Admixture K值范围 | Admixture K value range
    pca_components: int = 10  # PCA主成分数量 | PCA components count
    ld_window: int = 50  # LD窗口大小 | LD window size
    ld_step: int = 10  # LD步长 | LD step size
    ld_r2: float = 0.2  # LD r²阈值 | LD r² threshold

    # EMMAX分析参数 | EMMAX analysis parameters
    emmax_precision: int = 5  # 输出精度 | Output precision
    p_value_threshold: float = 5e-8  # 显著性阈值 | P-value threshold

    # 工作目录 | Working directory
    working_dir: str = "."

    # 内部属性 | Internal attributes
    base_name: str = "gwas_analysis"

    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.working_path = Path(self.working_dir).resolve()

        # 设置路径 | Setup paths
        self._setup_paths()

        # 标准化路径 | Normalize paths
        if self.vcf_file:
            self.vcf_file = os.path.normpath(os.path.abspath(self.vcf_file))

        if self.phenotype_file:
            self.phenotype_file = os.path.normpath(os.path.abspath(self.phenotype_file))

        self.manhattan_plot = os.path.normpath(os.path.abspath(self.manhattan_plot))
        self.qq_plot = os.path.normpath(os.path.abspath(self.qq_plot))
        self.significant_snps = os.path.normpath(os.path.abspath(self.significant_snps))
        self.excel_report = os.path.normpath(os.path.abspath(self.excel_report))
        self.html_report = os.path.normpath(os.path.abspath(self.html_report))

    def _setup_paths(self):
        """设置路径 | Setup paths"""
        # 如果没有指定输出文件，使用输出前缀 | If output files not specified, use output prefix
        if not self.manhattan_plot:
            self.manhattan_plot = f"{self.output_prefix}_manhattan.png"

        if not self.qq_plot:
            self.qq_plot = f"{self.output_prefix}_qq.png"

        if not self.significant_snps:
            self.significant_snps = f"{self.output_prefix}_significant_snps.txt"

        if not self.excel_report:
            self.excel_report = f"{self.output_prefix}_report.xlsx"

        if not self.html_report:
            self.html_report = f"{self.output_prefix}_report.html"

    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []

        # 检查输入文件 | Check input files
        if not self.vcf_file:
            errors.append("VCF文件路径必须指定 | VCF file path must be specified")
        elif not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在 | VCF file does not exist: {self.vcf_file}")

        if not self.phenotype_file:
            errors.append(
                "表型文件路径必须指定 | Phenotype file path must be specified"
            )
        elif not os.path.exists(self.phenotype_file):
            errors.append(
                f"表型文件不存在 | Phenotype file does not exist: {self.phenotype_file}"
            )

        # 检查输出目录 | Check output directories
        output_dirs = [
            (Path(self.manhattan_plot).parent, "曼哈顿图输出目录"),
            (Path(self.qq_plot).parent, "QQ图输出目录"),
            (Path(self.significant_snps).parent, "显著位点输出目录"),
            (Path(self.excel_report).parent, "Excel报告输出目录"),
            (Path(self.html_report).parent, "HTML报告输出目录"),
        ]

        for dir_path, desc in output_dirs:
            if not dir_path.exists():
                try:
                    dir_path.mkdir(parents=True, exist_ok=True)
                except Exception as e:
                    errors.append(
                        f"{desc}创建失败 | Failed to create {desc}: {dir_path}, 错误: {e}"
                    )

        # 检查参数范围 | Check parameter ranges
        if not 0 < self.maf_threshold < 0.5:
            errors.append(
                "MAF阈值必须在0-0.5之间 | MAF threshold must be between 0-0.5"
            )

        if not 0 < self.missing_threshold < 1:
            errors.append(
                "缺失率阈值必须在0-1之间 | Missing threshold must be between 0-1"
            )

        if not 0 < self.ld_r2 <= 1:
            errors.append(
                "LD r²阈值必须在0-1之间 | LD r² threshold must be between 0-1"
            )

        if self.pca_components < 1:
            errors.append(
                "PCA主成分数量必须大于0 | PCA components must be greater than 0"
            )

        if self.admixture_k_range[0] >= self.admixture_k_range[1]:
            errors.append("Admixture K值范围无效 | Invalid Admixture K value range")

        if errors:
            raise ValueError("\n".join(errors))

        return True
