"""
基因组组装质量评估配置管理模块|Genome Assembly Quality Control Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path, get_tool_path, resolve_legacy_path


@dataclass
class AssemblyQCConfig:
    """基因组组装质量评估配置类|Genome Assembly Quality Control Configuration Class"""

    # ==================== 必需输入参数 | Required Input Parameters ====================
    genome: str  # 基因组FASTA文件|Genome FASTA file
    lineage: str  # BUSCO谱系（如embryophyta_odb10）|BUSCO lineage (e.g., embryophyta_odb10)
    output_dir: str  # 输出目录|Output directory

    # ==================== 基因组统计参数 | Genome Statistics Parameters ====================
    sample_name: str = "genome_sample"  # 样品名称|Sample name
    threads: int = 12  # 全局线程数|Global threads (会自动分配给各子模块)

    # ==================== 核心评估：BUSCO | Core Evaluation: BUSCO ====================
    conda_env_busco: str = "~/miniforge3/envs/BUSCO_v.6.0.0"  # BUSCO conda环境路径|BUSCO conda env path
    busco_bin: str = field(
        default_factory=lambda: get_tool_path('busco', 'busco', 'BUSCO_BIN')
    )
    busco_dataset_path: str = "~/database/busco"  # BUSCO数据集路径|BUSCO dataset path
    busco_threads: int = 12  # BUSCO线程数|BUSCO threads
    busco_mode: str = "genome"  # BUSCO模式: genome/protein/transcriptome
    busco_download_path: Optional[str] = None  # BUSCO下载路径|BUSCO download path

    # ==================== 核心评估：LAI | Core Evaluation: LAI ====================
    lai_threads: int = 12  # LAI线程数|LAI threads
    lai_mode: str = "edta"  # LAI运行模式: edta/full/harvest/retrieve/calculate
    lai_quick_mode: bool = True  # LAI快速模式（使用-qq参数，跳过blastn）|LAI quick mode (use -qq flag, skip blastn) [default: True]

    # LAI工具路径|LAI tool paths
    conda_env_edta: str = "~/miniforge3/envs/EDTA_v.2.2.2"  # EDTA conda环境路径|EDTA conda env path
    conda_env_ltr_harvest: str = "~/miniforge3/envs/ltr_harvest_parallel_v.1.2"
    conda_env_ltr_finder: str = "~/miniforge3/envs/ltr_finder_parallel_v.1.3"
    conda_env_ltr_retriever: str = "~/miniforge3/envs/LTR_retriever_v.3.0.4"

    # ==================== 可选评估：QV | Optional Evaluation: QV ====================
    enable_qv: bool = True  # 启用QV评估|Enable QV evaluation [default: True]
    enable_qv_ngs: bool = True  # 启用NGS数据QV评估|Enable NGS QV evaluation [default: True]
    enable_qv_long_read: bool = True  # 启用long-read数据QV评估|Enable long-read QV evaluation [default: True]
    ngs_reads: Optional[str] = None  # NGS reads目录（用于QV和mapping）|NGS reads directory (for QV and mapping)
    long_reads: Optional[str] = None  # Long-reads目录（用于QV和mapping）|Long-reads directory (for QV and mapping)
    long_read_type: str = "hifi"  # Long-read数据类型: ont/pacbio/hifi|Long-read type
    qv_kmer_size: Optional[int] = None  # k-mer大小（None表示自动选择）|K-mer size (None for auto)
    qv_threads: int = 12  # QV线程数|QV threads

    # Merqury工具路径|Merqury tool paths
    conda_env_merqury: str = "~/miniforge3/envs/merqury_v.1.3"  # Merqury conda环境路径|Merqury conda env path
    merqury_path: Optional[str] = None  # Merqury安装路径（None表示从conda环境查找）
    qv_data_type: str = "auto"  # 数据类型: auto/illumina/hifi

    # ==================== 可选评估：Mapping | Optional Evaluation: Mapping ====================
    enable_mapping: bool = True  # 启用NGS Mapping评估|Enable NGS mapping evaluation [default: True]
    enable_long_read_mapping: bool = True  # 启用三代数据Mapping评估|Enable long-read mapping evaluation [default: True]
    mapping_threads: int = 12  # Mapping线程数|Mapping threads
    mapping_pattern: str = "_1.clean.fq.gz"  # FASTQ文件匹配模式（NGS）|FASTQ file pattern (NGS)
    long_read_mapping_pattern: str = "*.fq.gz"  # FASTQ文件匹配模式|FASTQ file pattern

    # Mapping工具conda环境路径|Mapping tool conda env paths
    conda_env_bwa: str = "~/miniforge3/envs/Population_genetics"  # BWA conda环境路径|BWA conda env path
    conda_env_samtools: str = "~/miniforge3/envs/GATK_v.4.6.2.0"  # Samtools conda环境路径|Samtools conda env path
    conda_env_minimap2: str = "~/miniforge3/envs/Genome_dedup"  # Minimap2 conda环境路径|Minimap2 conda env path

    # Mapping工具直接路径（优先使用，避免conda run的管道问题）|Mapping tool direct paths (preferred, avoid conda run pipe issues)
    bwa_path: str = "bwa"  # BWA可执行文件路径|BWA executable path (use system PATH or absolute path)
    samtools_path: str = "samtools"  # Samtools可执行文件路径|Samtools executable path
    minimap2_path: str = "minimap2"  # Minimap2可执行文件路径|Minimap2 executable path

    # Mapping工具路径|Mapping tool paths (deprecated, use conda_env_* instead)
    bwa_bin: str = field(
        default_factory=lambda: get_tool_path('bwa', 'bwa', 'BWA_BIN')
    )
    samtools_bin: str = field(
        default_factory=lambda: get_tool_path('samtools', 'samtools', 'SAMTOOLS_BIN')
    )
    bedtools_bin: str = field(
        default_factory=lambda: get_tool_path('bedtools', 'bedtools', 'BEDTOOLS_BIN')
    )

    # Mapping分析参数|Mapping analysis parameters
    min_mapping_quality: int = 20  # 最小比对质量|Minimum mapping quality
    min_base_quality: int = 20  # 最小碱基质量|Minimum base quality
    window_size: int = 1000000  # 滑窗大小（bp）|Window size for sliding window
    step_size: int = 100000  # 步长（bp）|Step size for sliding window

    # ==================== 优化：复用reads | Optimization: Reuse Reads ====================
    # 已移除reads_for_all，改用ngs_reads和long_reads统一管理
    # Removed reads_for_all, using ngs_reads and long_reads for unified management

    # ==================== 报告参数 | Report Parameters ====================
    generate_html: bool = True  # 生成HTML报告|Generate HTML report
    generate_table: bool = True  # 生成表格|Generate publication table
    table_format: str = "both"  # 表格格式: tsv/xlsx/both

    # ==================== 步骤控制 | Step Control ====================
    skip_busco: bool = False  # 跳过BUSCO评估|Skip BUSCO evaluation
    skip_lai: bool = False  # 跳过LAI评估|Skip LAI evaluation
    skip_qv: bool = False  # 跳过QV评估|Skip QV evaluation
    skip_mapping: bool = False  # 跳过Mapping评估|Skip mapping evaluation
    resume: bool = True  # 断点续传|Resume from previous run

    # ==================== 错误处理配置 | Error Handling Configuration ====================
    nfs_retry_max: int = 3  # NFS I/O错误最大重试次数|Max retries for NFS I/O errors
    nfs_retry_delay: int = 5  # NFS I/O错误重试间隔（秒）|Retry delay for NFS I/O errors (seconds)

    # ==================== 内部变量 | Internal Variables ====================
    work_dir: Optional[str] = None  # 工作目录|Work directory
    busco_output_dir: Optional[str] = None  # BUSCO输出目录|BUSCO output dir
    lai_output_dir: Optional[str] = None  # LAI输出目录|LAI output dir
    qv_output_dir: Optional[str] = None  # QV输出目录|QV output dir
    mapping_output_dir: Optional[str] = None  # Mapping输出目录|Mapping output dir

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有路径|Expand all paths
        # 关键：先展开~，再转换为绝对路径|CRITICAL: Expand ~ first, then convert to absolute path
        # 外部工具（BUSCO等）需要绝对路径|External tools (BUSCO, etc.) need absolute paths
        self.genome = os.path.abspath(expand_path(self.genome))
        self.output_dir = os.path.abspath(expand_path(self.output_dir))

        # 处理BUSCO lineage参数：如果是完整路径，提取数据集名称
        # Process BUSCO lineage parameter: extract dataset name if full path
        if os.path.isabs(self.lineage):
            # 从完整路径中提取数据集名称|Extract dataset name from full path
            self.lineage = os.path.basename(self.lineage)
        # lineage现在是数据集名称（如brassicales_odb12）
        # busco_dataset_path是数据集父目录
        # Now lineage is dataset name (e.g., brassicales_odb12)
        # busco_dataset_path is parent directory of datasets

        if self.ngs_reads:
            self.ngs_reads = os.path.abspath(expand_path(self.ngs_reads))

        if self.long_reads:
            self.long_reads = os.path.abspath(expand_path(self.long_reads))

        # 展开工具路径|Expand tool paths
        self.busco_bin = expand_path(self.busco_bin)
        self.bwa_bin = expand_path(self.bwa_bin)
        self.samtools_bin = expand_path(self.samtools_bin)
        self.bedtools_bin = expand_path(self.bedtools_bin)

        # 展开数据集路径|Expand dataset paths
        self.busco_dataset_path = expand_path(self.busco_dataset_path)

        # 展开conda环境路径|Expand conda environment paths
        self.conda_env_busco = expand_path(self.conda_env_busco)
        self.conda_env_edta = expand_path(self.conda_env_edta)
        self.conda_env_ltr_harvest = expand_path(self.conda_env_ltr_harvest)
        self.conda_env_ltr_finder = expand_path(self.conda_env_ltr_finder)
        self.conda_env_ltr_retriever = expand_path(self.conda_env_ltr_retriever)
        self.conda_env_merqury = expand_path(self.conda_env_merqury)
        self.conda_env_bwa = expand_path(self.conda_env_bwa)
        self.conda_env_samtools = expand_path(self.conda_env_samtools)
        self.conda_env_minimap2 = expand_path(self.conda_env_minimap2)

        if self.merqury_path:
            self.merqury_path = expand_path(self.merqury_path)

        # 创建输出目录|Create output directories
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        self.work_dir = self.output_dir

        # 定义各步骤输出目录|Define step output directories (遵循开发规范12.2)
        # 优先下划线规范名，回退点号老名用于断点续传|Prefer underscore, fall back to legacy dot name
        self.busco_output_dir = resolve_legacy_path(self.work_dir, "01_busco_evaluation")
        self.lai_output_dir = resolve_legacy_path(self.work_dir, "02_lai_evaluation")
        self.qv_output_dir = resolve_legacy_path(self.work_dir, "03_qv_evaluation")
        self.mapping_output_dir = resolve_legacy_path(self.work_dir, "04_mapping_evaluation")
        self.long_read_mapping_output_dir = resolve_legacy_path(self.work_dir, "05_long_read_mapping_evaluation")

        # 创建各步骤目录|Create step directories
        for step_dir in [
            self.busco_output_dir,
            self.lai_output_dir,
            self.qv_output_dir,
            self.mapping_output_dir,
            self.long_read_mapping_output_dir,
        ]:
            Path(step_dir).mkdir(parents=True, exist_ok=True)

        # 自动分配线程给各子模块|Automatically allocate threads to sub-modules
        self._allocate_threads()

        # 自动启用评估：如果提供了reads数据，自动启用对应的评估
        # Auto-enable evaluations: if reads data provided, auto-enable corresponding evaluations
        if self.ngs_reads:
            self.enable_qv_ngs = True
            self.enable_mapping = True

        if self.long_reads:
            self.enable_qv_long_read = True
            self.enable_long_read_mapping = True

        # 总体QV开关：只要有任何一种reads数据，就启用QV
        # Overall QV switch: enable QV if any reads data provided
        if self.ngs_reads or self.long_reads:
            self.enable_qv = True

    def _allocate_threads(self):
        """自动分配线程给各子模块|Automatically allocate threads to sub-modules"""
        # 由于工作流是串行执行的，每个步骤都使用全部线程以最大化性能
        # 如果用户明确指定了各子模块的线程数，则使用用户指定的值

        # BUSCO: 使用全部线程
        if not hasattr(self, '_user_specified_busco_threads'):
            self.busco_threads = max(1, self.threads)

        # LAI: 使用全部线程（串行执行，不需要分配）
        if not hasattr(self, '_user_specified_lai_threads'):
            self.lai_threads = max(1, self.threads)

        # QV: 使用全部线程（串行执行，不需要分配）
        if not hasattr(self, '_user_specified_qv_threads'):
            self.qv_threads = max(1, self.threads)

        # Mapping: 使用全部线程
        if not hasattr(self, '_user_specified_mapping_threads'):
            self.mapping_threads = max(1, self.threads)

        # 三代数据Mapping: 使用全部线程
        if not hasattr(self, '_user_specified_long_read_mapping_threads'):
            self.long_read_mapping_threads = max(1, self.threads)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome}")

        # 检查BUSCO数据集路径|Check BUSCO dataset path
        if not os.path.exists(self.busco_dataset_path):
            errors.append(f"BUSCO数据集路径不存在|BUSCO dataset path not found: {self.busco_dataset_path}")

        # 检查QV评估参数|Check QV evaluation parameters
        if self.enable_qv:
            if not self.ngs_reads and not self.long_reads:
                errors.append("启用QV评估时需要提供ngs_reads或long_reads|ngs_reads or long_reads required when QV enabled")

        # 检查Mapping评估参数|Check mapping evaluation parameters
        if self.enable_mapping and not self.ngs_reads:
            errors.append("启用NGS Mapping评估时需要提供ngs_reads|ngs_reads required when NGS mapping enabled")

        # 检查三代数据Mapping评估参数|Check long-read mapping evaluation parameters
        if self.enable_long_read_mapping and not self.long_reads:
            errors.append("启用三代数据Mapping评估时需要提供long_reads|long_reads required when long-read mapping enabled")

        # 检查参数范围|Check parameter ranges
        if self.busco_threads <= 0:
            errors.append(f"BUSCO线程数必须大于0|BUSCO threads must be > 0: {self.busco_threads}")

        if self.lai_threads <= 0:
            errors.append(f"LAI线程数必须大于0|LAI threads must be > 0: {self.lai_threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_mapping_reads(self) -> Optional[str]:
        """获取用于Mapping评估的reads目录|Get reads directory for mapping evaluation"""
        return self.ngs_reads

    def get_long_read_mapping_reads(self) -> Optional[str]:
        """获取用于三代数据Mapping评估的reads目录|Get reads directory for long-read mapping evaluation"""
        return self.long_reads

    def get_qv_ngs_reads(self) -> Optional[str]:
        """获取用于NGS QV评估的reads目录|Get reads directory for NGS QV evaluation"""
        return self.ngs_reads

    def get_qv_long_read_reads(self) -> Optional[str]:
        """获取用于long-read QV评估的reads目录|Get reads directory for long-read QV evaluation"""
        return self.long_reads

    def get_summary(self) -> dict:
        """获取配置摘要|Get configuration summary"""
        return {
            # 样本信息|Sample info
            "sample_name": self.sample_name,
            "genome": self.genome,
            "lineage": self.lineage,

            # 核心评估|Core evaluation
            "busco_enabled": not self.skip_busco,
            "lai_enabled": not self.skip_lai,

            # 可选评估|Optional evaluation
            "qv_enabled": self.enable_qv and not self.skip_qv,
            "mapping_enabled": self.enable_mapping and not self.skip_mapping,

            # 输出目录|Output directories
            "output_dir": self.output_dir,
            "busco_output_dir": self.busco_output_dir,
            "lai_output_dir": self.lai_output_dir,
            "qv_output_dir": self.qv_output_dir,
            "mapping_output_dir": self.mapping_output_dir,
        }
