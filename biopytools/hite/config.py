"""
HiTE 配置管理模块|HiTE Configuration Management Module

包含 HiTE 和 panHiTE 的配置数据类
Contains configuration dataclasses for HiTE and panHiTE
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import os
from ..common.paths import expand_path


@dataclass
class HiteConfig:
    """
    HiTE 单基因组配置类|HiTE Single-genome Configuration Class

    用于 HiTE 单基因组转座子检测与注释的配置管理
    Configuration management for HiTE single-genome transposon detection and annotation
    """

    # ==================== 必需参数|Required parameters ====================
    genome: str
    """基因组FASTA文件路径|Genome FASTA file path"""

    # ==================== 路径配置|Path configuration ====================
    singularity_cmd: str = '~/miniforge3/envs/singularity_v.3.8.7/bin/singularity'
    """Singularity可执行文件路径|Singularity executable path"""

    sif_file: str = '~/software/singularity/hite_3.3.3.sif'
    """HiTE SIF镜像文件路径|HiTE SIF image file path"""

    output_dir: str = './hite_output'
    """输出目录路径|Output directory path"""

    work_dir: str = '/tmp'
    """临时工作目录路径|Temporary work directory path"""

    # ==================== 处理参数|Processing parameters ====================
    threads: int = 12
    """线程数|Number of threads"""

    plant: bool = True
    """是否为植物基因组|Whether the genome is from plant"""

    annotate: bool = False
    """是否注释基因组|Whether to annotate the genome with detected TEs"""

    recover: bool = False
    """是否启用断点续跑|Whether to enable recovery mode"""

    domain: bool = False
    """是否预测TE蛋白结构域|Whether to predict protein domains in TEs"""

    # ==================== 高级参数|Advanced parameters ====================
    chunk_size: int = 400
    """基因组分块大小(MB)|Genome chunk size in MB"""

    miu: float = 1.3e-8
    """中性突变率(per bp per year)|Neutral mutation rate"""

    min_te_len: int = 80
    """最小TE长度(bp)|Minimum TE length in bp"""

    te_type: str = 'all'
    """检测的TE类型|TE type to detect: [ltr|tir|helitron|non-ltr|all]"""

    remove_nested: bool = True
    """是否移除嵌套TE|Whether to remove nested TEs"""

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
        self.genome = os.path.abspath(self.genome)
        self.output_dir = os.path.abspath(self.output_dir)
        self.work_dir = os.path.abspath(self.work_dir)

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

        # 检查必需文件
        if not os.path.exists(self.genome):
            errors.append(
                f"基因组文件不存在|Genome file not found: {self.genome}"
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

        if self.min_te_len < 1:
            errors.append(
                f"最小TE长度必须大于0|Min TE length must be greater than 0: "
                f"{self.min_te_len}"
            )

        if self.chunk_size < 1:
            errors.append(
                f"分块大小必须大于0|Chunk size must be greater than 0: "
                f"{self.chunk_size}"
            )

        if self.miu <= 0:
            errors.append(
                f"突变率必须大于0|Mutation rate must be greater than 0: "
                f"{self.miu}"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return True


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
