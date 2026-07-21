"""
HiTE 配置管理模块|HiTE Configuration Management Module

包含 HiTE 和 panHiTE 的配置数据类
Contains configuration dataclasses for HiTE and panHiTE
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import os
from ..common.paths import get_tool_path, expand_path


@dataclass
class HiteConfig:
    """
    HiTE 单基因组配置类|HiTE Single-genome Configuration Class

    singularity 直接挂载模式:singularity exec --bind host:host <sif> python /HiTE/main.py
    Singularity direct-mount mode: HiTE writes directly to host output_dir/01_hite
    """

    # ==================== 必需参数|Required ====================
    genome: str
    """基因组FASTA文件路径|Genome FASTA file path"""

    # ==================== 输出|Output ====================
    output_dir: str = "./hite_output"
    """输出根目录|Output root directory"""

    # ==================== singularity 容器|Container ====================
    singularity_path: str = field(
        default_factory=lambda: get_tool_path(
            "singularity",
            "~/miniforge3/envs/singularity_v.3.8.7/bin/singularity",
            "SINGULARITY_PATH",
        )
    )
    """Singularity 可执行文件|Singularity executable"""

    sif_file: str = field(
        default_factory=lambda: get_tool_path(
            "hite",
            "~/software/singularity/hite_3.3.3.sif",
            "HITE_SIF",
        )
    )
    """HiTE SIF 镜像|HiTE SIF image"""

    # ==================== 处理参数|Processing ====================
    threads: int = 12
    plant: bool = True
    annotate: bool = False
    recover: bool = False
    domain: bool = False

    # ==================== 高级参数|Advanced ====================
    chunk_size: int = 400
    miu: float = 1.3e-8
    min_te_len: int = 80
    te_type: str = "all"
    remove_nested: bool = True
    curated_lib: Optional[str] = None
    debug: bool = False

    def __post_init__(self):
        """展开所有 ~ 路径并创建目录|Expand all ~ paths and create dirs"""
        # 关键:所有含 ~ 路径必须在此展开|CRITICAL: expand all ~ paths here
        self.singularity_path = expand_path(self.singularity_path)
        self.sif_file = expand_path(self.sif_file)
        self.genome = expand_path(self.genome)
        self.output_dir = expand_path(self.output_dir)
        if self.curated_lib:
            self.curated_lib = expand_path(self.curated_lib)

        # 派生目录|Derived directories (规范 §12 by-step 结构)
        self.hite_out_dir = os.path.join(self.output_dir, "01_hite")
        self.work_dir = os.path.join(self.hite_out_dir, "work")
        self.pipeline_info_dir = os.path.join(self.output_dir, "00_pipeline_info")
        self.logs_dir = os.path.join(self.output_dir, "99_logs")

        # 创建目录|Create directories
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        Path(self.hite_out_dir).mkdir(parents=True, exist_ok=True)
        Path(self.pipeline_info_dir).mkdir(parents=True, exist_ok=True)
        Path(self.logs_dir).mkdir(parents=True, exist_ok=True)

    def validate(self) -> bool:
        """验证配置|Validate configuration"""
        errors = []

        if not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome}")
        if not os.path.exists(self.singularity_cmd_path_check()):
            errors.append(
                f"Singularity不存在|Singularity not found: {self.singularity_path}"
            )
        if not os.path.exists(self.sif_file):
            errors.append(f"HiTE SIF文件不存在|HiTE SIF file not found: {self.sif_file}")

        if self.te_type not in ['ltr', 'tir', 'helitron', 'non-ltr', 'all']:
            errors.append(
                f"无效的TE类型|Invalid TE type: {self.te_type} "
                f"(应为|should be ltr|tir|helitron|non-ltr|all)"
            )
        if self.threads < 1:
            errors.append(f"线程数必须大于0|Threads must be > 0: {self.threads}")
        if self.min_te_len < 1:
            errors.append(f"最小TE长度必须大于0|Min TE length must be > 0: {self.min_te_len}")
        if self.chunk_size < 1:
            errors.append(f"分块大小必须大于0|Chunk size must be > 0: {self.chunk_size}")
        if self.miu <= 0:
            errors.append(f"突变率必须大于0|Mutation rate must be > 0: {self.miu}")

        if errors:
            raise ValueError("\n".join(errors))
        return True

    def singularity_cmd_path_check(self) -> str:
        """返回 singularity 路径供存在性检查|Return singularity path for existence check"""
        return self.singularity_path


@dataclass
class PanHiteConfig:
    """
    panHiTE 群体基因组配置类|panHiTE Pan-genome Configuration Class

    用于 panHiTE 泛基因组转座子检测与比较分析的配置管理
    Configuration management for panHiTE pan-genome transposon detection
    and comparative analysis
    """

    # ==================== 必需参数|Required parameters ====================
    pan_genomes_dir: str
    """泛基因组目录路径|Pan-genomes directory path"""

    genome_list: str
    """基因组列表文件路径|Genome list file path"""

    # ==================== 可选路径|Optional paths ====================
    genes_dir: Optional[str] = None
    """基因注释文件目录路径|Gene annotation files directory path"""

    rna_dir: Optional[str] = None
    """RNA-seq数据目录路径|RNA-seq data directory path"""

    # ==================== 路径配置|Path configuration ====================
    singularity_cmd: str = '~/miniforge3/envs/singularity_v.3.8.7/bin/singularity'
    """Singularity可执行文件路径|Singularity executable path"""

    sif_file: str = '~/software/singularity/hite_3.3.3.sif'
    """HiTE SIF镜像文件路径|HiTE SIF image file path"""

    output_dir: str = './panhite_output'
    """输出目录路径|Output directory path"""

    # ==================== 处理参数|Processing parameters ====================
    threads: int = 12
    """线程数|Number of threads"""

    miu: float = 1.3e-8
    """中性突变率(per bp per year)|Neutral mutation rate"""

    te_type: str = 'all'
    """检测的TE类型|TE type to detect: [ltr|tir|helitron|non-ltr|all]"""

    skip_analyze: bool = False
    """是否跳过分析，仅生成panTE库|Whether to skip analysis, only generate panTE library"""

    recover: bool = False
    """是否启用断点续跑|Whether to enable recovery mode"""

    debug: bool = False
    """是否开启调试模式|Whether to enable debug mode"""

    # ==================== 容器配置|Container configuration ====================
    mount_home: bool = True
    """是否挂载用户主目录|Whether to mount user home directory"""

    def __post_init__(self):
        """
        初始化后处理|Post-initialization processing

        标准化所有路径并创建输出目录
        Normalize all paths and create output directory
        """
        # 展开 ~ 为完整路径
        self.singularity_cmd = os.path.expanduser(self.singularity_cmd)
        self.sif_file = os.path.abspath(self.sif_file)
        self.pan_genomes_dir = os.path.abspath(self.pan_genomes_dir)
        self.genome_list = os.path.abspath(self.genome_list)
        self.output_dir = os.path.abspath(self.output_dir)

        # 处理可选路径
        if self.genes_dir:
            self.genes_dir = os.path.abspath(self.genes_dir)
        if self.rna_dir:
            self.rna_dir = os.path.abspath(self.rna_dir)

        # 创建输出目录
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """
        验证配置参数|Validate configuration parameters

        Returns:
            bool: 验证成功返回True|Returns True if validation succeeds

        Raises:
            ValueError: 如果验证失败|If validation fails
        """
        errors = []

        # 检查必需目录
        if not os.path.exists(self.pan_genomes_dir):
            errors.append(
                f"泛基因组目录不存在|Pan-genomes directory not found: "
                f"{self.pan_genomes_dir}"
            )

        if not os.path.exists(self.genome_list):
            errors.append(
                f"基因组列表文件不存在|Genome list file not found: "
                f"{self.genome_list}"
            )

        # 检查可选目录
        if self.genes_dir and not os.path.exists(self.genes_dir):
            errors.append(
                f"基因注释目录不存在|Genes directory not found: {self.genes_dir}"
            )

        if self.rna_dir and not os.path.exists(self.rna_dir):
            errors.append(
                f"RNA-seq目录不存在|RNA-seq directory not found: {self.rna_dir}"
            )

        # 检查Singularity
        if not os.path.exists(self.singularity_cmd):
            errors.append(
                f"Singularity不存在|Singularity not found: {self.singularity_cmd}"
            )

        # 检查SIF文件
        if not os.path.exists(self.sif_file):
            errors.append(
                f"HiTE SIF文件不存在|HiTE SIF file not found: {self.sif_file}"
            )

        # 验证参数范围
        if self.te_type not in ['ltr', 'tir', 'helitron', 'non-ltr', 'all']:
            errors.append(
                f"无效的TE类型|Invalid TE type: {self.te_type} "
                f"(应为|should be ltr|tir|helitron|non-ltr|all)"
            )

        if self.threads < 1:
            errors.append(
                f"线程数必须大于0|Threads must be greater than 0: {self.threads}"
            )

        if self.miu <= 0:
            errors.append(
                f"突变率必须大于0|Mutation rate must be greater than 0: {self.miu}"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return True
