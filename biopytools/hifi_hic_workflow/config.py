"""
HiFi+Hi-C工作流配置管理模块|HiFi+Hi-C Workflow Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path, get_tool_path


@dataclass
class HifiHicWorkflowConfig:
    """HiFi+Hi-C工作流配置类|HiFi+Hi-C Workflow Configuration Class"""

    # ==================== 必需输入参数 | Required Input Parameters ====================
    hifi_reads: str  # HiFi reads文件|HiFi reads file
    hic_r1: str  # Hi-C R1文件|Hi-C R1 file
    hic_r2: str  # Hi-C R2文件|Hi-C R2 file
    reference_genome: str  # 参考基因组（仅用于命名）|Reference genome (for naming only)
    output_dir: str  # 输出目录|Output directory

    # ==================== 全局参数 | Global Parameters ====================
    prefix: str = "genome_sample"  # 样本前缀|Sample prefix
    threads: int = 64  # 默认线程数|Default threads

    # ==================== 流程控制 | Workflow Control ====================
    skip_hifi_hic: bool = False  # 跳过HiFi组装|Skip HiFi assembly
    skip_haphic: bool = False  # 跳过Hi-C挂载|Skip Hi-C scaffolding
    skip_rename: bool = False  # 跳过重命名|Skip renaming
    skip_heatmap: bool = False  # 跳过热图|Skip heatmap
    resume: bool = True  # 断点续传|Resume from previous run
    force_rerun: bool = False  # 强制重新运行|Force rerun all steps

    # ==================== Step 1: HiFi组装参数 | Step 1: HiFi Assembly Parameters ====================
    # 基本参数|Basic parameters
    genome_size: str = "1.45g"  # 基因组大小|Genome size estimate
    n_hap: int = 2  # 倍性|Ploidy
    purge_level: Optional[int] = None  # purge level (-l): 0=no purging, 1=light, 2/3=aggressive
    hom_cov: Optional[int] = None  # homozygous read coverage (--hom-cov)

    # NGS Polish参数|NGS polish parameters
    use_ngs_polish: bool = False  # 是否使用NGS polish|Whether to use NGS polish
    ngs_data: Optional[str] = None  # NGS二代数据目录|NGS second-generation data directory
    ngs_high_cov: float = 95.0  # 高质量contig覆盖度阈值|High quality contig coverage threshold
    ngs_pattern: str = "_1.clean.fq.gz"  # NGS文件匹配模式|NGS file matching pattern

    # ==================== Step 2: HapHiC挂载参数 | Step 2: HapHiC Scaffolding Parameters ====================
    nchrs: Optional[int] = None  # 染色体数量（如果为None，从reference_genome统计）|Number of chromosomes (count from reference if None)

    # HapHiC工具路径|HapHiC tool paths
    haphic_bin: str = field(
        default_factory=lambda: get_tool_path('haphic', 'haphic', 'HAPHIC_BIN')
    )
    bwa_bin: str = field(
        default_factory=lambda: get_tool_path('bwa', 'bwa', 'BWA_BIN')
    )
    samtools_bin: str = field(
        default_factory=lambda: get_tool_path('samtools', 'samtools', 'SAMTOOLS_BIN')
    )

    # Hi-C数据处理参数|Hi-C data processing parameters
    mapq_threshold: int = 1  # MAPQ阈值|MAPQ threshold
    edit_distance: int = 3  # 编辑距离|Edit distance
    min_re_sites: int = 25  # 最小限制性酶切位点数|Minimum restriction sites

    # 聚类参数|Clustering parameters
    min_inflation: float = 1.0  # 最小膨胀参数|Min inflation
    max_inflation: float = 3.0  # 最大膨胀参数|Max inflation
    inflation_step: float = 0.2  # 膨胀参数步长|Inflation step
    nx: int = 80  # 聚类参数|Clustering parameter
    min_group_len: int = 0  # 最小组长度|Minimum group length

    # 排序和定向参数|Ordering and orientation parameters
    processes: int = 8  # 进程数|Number of processes
    fast_sorting: bool = True  # 快速排序|Fast sorting
    allhic_optimization: bool = False  # ALLHiC优化|ALLHiC optimization

    # 组装校正参数|Assembly correction parameters
    correct_nrounds: int = 2  # 校正轮数|Correction rounds
    correct_min_coverage: float = 10.0  # 校正最小覆盖度|Correction min coverage
    correct_resolution: Optional[int] = None  # 校正分辨率|Correction resolution

    # 可视化参数|Visualization parameters
    generate_plots: bool = False  # 生成图表|Generate plots
    bin_size: int = 500  # 装箱大小|Bin size
    min_len: float = 1.0  # 最小长度|Minimum length
    separate_plots: bool = False  # 分离图表|Separate plots

    # 性能参数|Performance parameters
    haphic_threads: int = 64  # HapHiC线程数|HapHiC threads
    memory_limit: Optional[str] = "300G"  # 内存限制|Memory limit

    # 高级选项|Advanced options
    quick_view: bool = False  # 快速查看|Quick view
    re_sites: str = "GATC"  # 限制性内切酶位点|Restriction enzyme sites
    skip_clustering: bool = False  # 跳过聚类|Skip clustering

    # ==================== Step 3: 染色体重命名参数 | Step 3: Chromosome Rename Parameters ====================
    rename_keep_all: bool = True  # 保留所有序列|Keep all sequences (chrNN + scaffolds)

    # 参考基因组命名参数|Reference genome naming parameters
    naming_min_identity: float = 80.0  # 最小序列一致性|Min sequence identity (%)
    naming_min_coverage: float = 80.0  # 最小覆盖度|Min coverage (%)
    naming_minimap2_preset: str = "asm5"  # minimap2预设|minimap2 preset (asm5/asm10/asm20)

    # ==================== Step 4: Hi-C热图参数 | Step 4: Hi-C Heatmap Parameters ====================
    # HiCPro参数|HiCPro parameters
    hicpro_restriction_enzyme: str = "MboI"  # 限制性内切酶|Restriction enzyme
    hicpro_bin_sizes: str = "20000 40000 150000 500000 1000000"  # bin大小列表|Bin sizes
    hicpro_max_memory_gb: int = 200  # HiCPro最大内存|HiCPro max memory
    hicpro_matrix_format: str = "upper"  # 矩阵格式|Matrix format (upper/lower/complete)

    # PlotHiC参数|PlotHiC parameters
    heatmap_resolution: int = 100000  # 热图分辨率|Heatmap resolution (bp)
    heatmap_color_map: str = "YlOrRd"  # 颜色方案|Color scheme
    heatmap_dpi: int = 300  # 图像分辨率|Image resolution
    heatmap_format: str = "pdf"  # 输出格式|Output format
    heatmap_bar_max: int = 1  # 颜色条最大值|Color bar max value

    # Hi-C热图工具路径|Hi-C heatmap tool paths
    hicpro_path: str = field(
        default_factory=lambda: get_tool_path(
            'hicpro',
            '~/software/HiC-Pro_v3.1.0/HiC-Pro_3.1.0/bin/HiC-Pro',
            'HICPRO_PATH'
        )
    )
    hicpro_sif: Optional[str] = field(
        default_factory=lambda: get_tool_path(
            'hicpro_sif',
            '',
            'HICPRO_SIF'
        )
    )
    plothic_path: str = field(
        default_factory=lambda: get_tool_path(
            'plothic',
            '~/miniforge3/envs/plothic_v.1.0.0/bin/plothic',
            'PLOTHIC_PATH'
        )
    )

    # ==================== 内部变量 | Internal Variables ====================
    work_dir: Optional[str] = None  # 工作目录|Work directory
    hifi_hic_output_dir: Optional[str] = None  # HiFi组装输出目录|HiFi assembly output dir
    haphic_output_dir: Optional[str] = None  # HapHiC输出目录|HapHiC output dir
    rename_output_dir: Optional[str] = None  # 重命名输出目录|Rename output dir
    heatmap_output_dir: Optional[str] = None  # 热图输出目录|Heatmap output dir

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有路径|Expand all paths
        self.hifi_reads = expand_path(self.hifi_reads)
        self.hic_r1 = expand_path(self.hic_r1)
        self.hic_r2 = expand_path(self.hic_r2)
        self.reference_genome = expand_path(self.reference_genome)
        self.output_dir = expand_path(self.output_dir)

        if self.ngs_data:
            self.ngs_data = expand_path(self.ngs_data)

        # 展开工具路径|Expand tool paths
        self.haphic_bin = expand_path(self.haphic_bin)
        self.bwa_bin = expand_path(self.bwa_bin)
        self.samtools_bin = expand_path(self.samtools_bin)
        self.hicpro_path = expand_path(self.hicpro_path)
        self.plothic_path = expand_path(self.plothic_path)

        if self.hicpro_sif:
            self.hicpro_sif = expand_path(self.hicpro_sif)

        # 创建工作目录|Create work directory
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        self.work_dir = self.output_dir

        # 定义各步骤输出目录|Define step output directories (遵循开发规范12.2)
        self.hifi_hic_output_dir = os.path.join(self.work_dir, "01.hifi_assembly")
        self.haphic_output_dir = os.path.join(self.work_dir, "02.hic_scaffolding")
        self.rename_output_dir = os.path.join(self.work_dir, "03.chromosome_rename")
        self.heatmap_output_dir = os.path.join(self.work_dir, "04.hic_heatmap")

        # 创建各步骤目录|Create step directories
        for step_dir in [
            self.hifi_hic_output_dir,
            self.haphic_output_dir,
            self.rename_output_dir,
            self.heatmap_output_dir
        ]:
            Path(step_dir).mkdir(parents=True, exist_ok=True)

        # 如果未指定nchrs，从reference_genome统计|Count nchrs from reference if not specified
        if self.nchrs is None:
            self.nchrs = self._count_reference_chromosomes()

    def _count_reference_chromosomes(self) -> int:
        """
        从参考基因组统计染色体数量|Count chromosomes from reference genome

        Returns:
            int: 染色体数量|Number of chromosomes
        """
        try:
            count = 0
            with open(self.reference_genome, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        count += 1
            return count
        except Exception as e:
            raise ValueError(f"无法读取参考基因组|Cannot read reference genome: {e}")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input files
        if not os.path.exists(self.hifi_reads):
            errors.append(f"HiFi reads文件不存在|HiFi reads file not found: {self.hifi_reads}")

        if not os.path.exists(self.hic_r1):
            errors.append(f"Hi-C R1文件不存在|Hi-C R1 file not found: {self.hic_r1}")

        if not os.path.exists(self.hic_r2):
            errors.append(f"Hi-C R2文件不存在|Hi-C R2 file not found: {self.hic_r2}")

        if not os.path.exists(self.reference_genome):
            errors.append(f"参考基因组文件不存在|Reference genome not found: {self.reference_genome}")

        # 检查NGS数据|Check NGS data
        if self.use_ngs_polish and self.ngs_data:
            if not os.path.exists(self.ngs_data):
                errors.append(f"NGS数据目录不存在|NGS data directory not found: {self.ngs_data}")
            if not os.path.isdir(self.ngs_data):
                errors.append(f"NGS数据路径不是目录|NGS data path is not a directory: {self.ngs_data}")

        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须大于0|Threads must be > 0: {self.threads}")

        if self.n_hap <= 0:
            errors.append(f"倍性必须大于0|Ploidy must be > 0: {self.n_hap}")

        if self.nchrs <= 0:
            errors.append(f"染色体数必须大于0|Chromosomes must be > 0: {self.nchrs}")

        if self.naming_min_identity <= 0 or self.naming_min_identity > 100:
            errors.append(f"命名一致性必须在0-100之间|Naming identity must be 0-100: {self.naming_min_identity}")

        if self.naming_min_coverage <= 0 or self.naming_min_coverage > 100:
            errors.append(f"命名覆盖度必须在0-100之间|Naming coverage must be 0-100: {self.naming_min_coverage}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_summary(self) -> dict:
        """
        获取配置摘要|Get configuration summary

        Returns:
            dict: 配置摘要|Configuration summary
        """
        return {
            # 输入文件|Input files
            "hifi_reads": self.hifi_reads,
            "hic_r1": self.hic_r1,
            "hic_r2": self.hic_r2,
            "reference_genome": self.reference_genome,
            "ngs_data": self.ngs_data if self.use_ngs_polish else None,
            # 全局参数|Global parameters
            "prefix": self.prefix,
            "threads": self.threads,
            "nchrs": self.nchrs,
            # 流程控制|Workflow control
            "skip_hifi_hic": self.skip_hifi_hic,
            "skip_haphic": self.skip_haphic,
            "skip_rename": self.skip_rename,
            "skip_heatmap": self.skip_heatmap,
            "use_ngs_polish": self.use_ngs_polish,
            "resume": self.resume,
            "force_rerun": self.force_rerun,
            # 输出目录|Output directories
            "work_dir": self.work_dir,
            "hifi_hic_output_dir": self.hifi_hic_output_dir,
            "haphic_output_dir": self.haphic_output_dir,
            "rename_output_dir": self.rename_output_dir,
            "heatmap_output_dir": self.heatmap_output_dir,
        }
