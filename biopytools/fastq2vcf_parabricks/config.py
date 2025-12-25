"""
Fastq到VCF (Parabricks) 配置管理模块 | Fastq to VCF (Parabricks) Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class Fastq2VcfParabricksConfig:
    """Fastq到VCF (Parabricks) 配置类 | Fastq to VCF (Parabricks) Configuration Class"""

    # 必需输入 | Required inputs
    raw_fastq_dir: str
    ref_genome_fa: str
    project_base: str

    # 路径配置 | Path configuration
    clean_fastq_dir: Optional[str] = None
    mapping_dir: Optional[str] = None
    gvcf_dir: Optional[str] = None
    bam_dir: Optional[str] = None
    joint_dir: Optional[str] = None
    filter_dir: Optional[str] = None
    output_dir: Optional[str] = None

    # 线程配置 | Thread configuration
    threads_mapping: int = 88
    threads_gtx: int = 88
    threads_filter: int = 88

    # 过滤参数 | Filtering parameters
    snp_min_dp: int = 5
    snp_min_qual: int = 30
    indel_min_dp: int = 5
    indel_min_qual: int = 30

    # 样本阈值 | Sample thresholds
    gatk_threshold: int = 4        # < 4 使用 GATK
    gtx_single_threshold: int = 200 # < 200 使用 GTX 单机模式
    gtx_window_size: int = 20000000 # GTX 分块窗口大小 (20Mb)

    # GTX WGS参数 | GTX WGS parameters
    gtx_pcr_indel_model: str = "CONSERVATIVE"
    gtx_min_confidence: int = 30
    gtx_min_base_qual: int = 20
    gtx_ploidy: int = 2

    # 工具路径 | Tool paths
    gtx_bin: str = "/share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx"

    # 高级选项 | Advanced options
    enable_checkpoint: bool = True
    dry_run: bool = False
    verbose: bool = False
    skip_qc: bool = False
    skip_mapping: bool = False
    read1_pattern_fastp: str = "_1.fq.gz"
    read2_pattern_fastp: str = "_2.fq.gz"

    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""

        # 设置默认路径 | Set default paths
        if self.clean_fastq_dir is None:
            self.clean_fastq_dir = os.path.join(self.project_base, "01.data", "clean")
        if self.mapping_dir is None:
            self.mapping_dir = os.path.join(self.project_base, "02.mapping")
        if self.gvcf_dir is None:
            self.gvcf_dir = os.path.join(self.mapping_dir, "vcf")
        if self.bam_dir is None:
            self.bam_dir = os.path.join(self.mapping_dir, "bam")
        if self.joint_dir is None:
            self.joint_dir = os.path.join(self.project_base, "03.joint_calling")
        if self.filter_dir is None:
            self.filter_dir = os.path.join(self.project_base, "04.filtered_snp_indel")
        if self.output_dir is None:
            self.output_dir = self.filter_dir

        # 创建输出路径对象 | Create output path objects
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径 | Normalize paths
        self.raw_fastq_dir = os.path.normpath(os.path.abspath(self.raw_fastq_dir))
        self.ref_genome_fa = os.path.normpath(os.path.abspath(self.ref_genome_fa))
        self.project_base = os.path.normpath(os.path.abspath(self.project_base))
        self.clean_fastq_dir = os.path.normpath(os.path.abspath(self.clean_fastq_dir))
        self.mapping_dir = os.path.normpath(os.path.abspath(self.mapping_dir))
        self.gvcf_dir = os.path.normpath(os.path.abspath(self.gvcf_dir))
        self.bam_dir = os.path.normpath(os.path.abspath(self.bam_dir))
        self.joint_dir = os.path.normpath(os.path.abspath(self.joint_dir))
        self.filter_dir = os.path.normpath(os.path.abspath(self.filter_dir))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        self.gtx_bin = os.path.normpath(os.path.abspath(self.gtx_bin))

        # 创建必要的目录 | Create necessary directories
        for dir_path in [self.clean_fastq_dir, self.mapping_dir, self.gvcf_dir,
                        self.bam_dir, self.joint_dir, self.filter_dir]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []

        # 检查必需的输入路径 | Check required input paths
        if not os.path.exists(self.raw_fastq_dir):
            errors.append(f"原始FASTQ目录不存在 | Raw FASTQ directory does not exist: {self.raw_fastq_dir}")

        if not os.path.exists(self.ref_genome_fa):
            errors.append(f"参考基因组文件不存在 | Reference genome file does not exist: {self.ref_genome_fa}")

        # 检查目录是否为空 | Check if directories are empty
        if os.path.exists(self.raw_fastq_dir):
            if not any(os.scandir(self.raw_fastq_dir)):
                errors.append(f"原始FASTQ目录为空 | Raw FASTQ directory is empty: {self.raw_fastq_dir}")

        # 检查GTX工具路径（如果需要） | Check GTX tool path (if needed)
        if not os.path.exists(self.gtx_bin):
            errors.append(f"GTX工具不存在 | GTX tool does not exist: {self.gtx_bin}")

        # 检查线程数是否合理 | Check if thread counts are reasonable
        if self.threads_mapping < 1:
            errors.append("比对线程数必须大于0 | Mapping threads must be greater than 0")
        if self.threads_gtx < 1:
            errors.append("GTX线程数必须大于0 | GTX threads must be greater than 0")
        if self.threads_filter < 1:
            errors.append("过滤线程数必须大于0 | Filtering threads must be greater than 0")

        # 检查过滤参数 | Check filtering parameters
        if self.snp_min_dp < 1:
            errors.append("SNP最小深度必须大于0 | SNP minimum depth must be greater than 0")
        if self.snp_min_qual < 0:
            errors.append("SNP最小质量不能为负数 | SNP minimum quality cannot be negative")
        if self.indel_min_dp < 1:
            errors.append("InDel最小深度必须大于0 | InDel minimum depth must be greater than 0")
        if self.indel_min_qual < 0:
            errors.append("InDel最小质量不能为负数 | InDel minimum quality cannot be negative")

        if errors:
            raise ValueError("\n".join(errors))

        return True