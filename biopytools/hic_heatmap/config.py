"""Hi-C热图分析配置模块|Hi-C heatmap analysis configuration module"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import os
from ..common.paths import expand_path, get_tool_path


@dataclass
class HiCConfig:
    """Hi-C热图配置类|Hi-C heatmap configuration class (Juicer + PlotHiC)"""

    # 必需参数|Required parameters
    genome: str
    output_dir: str
    genome_id: Optional[str] = None  # 基因组ID|Genome ID (for Juicer, auto-extracted if not provided)

    # 输入文件（三种模式）|Input files (three modes)
    fastq_r1: Optional[str] = None  # 模式1：从FASTQ开始|Mode 1: Start from FASTQ
    fastq_r2: Optional[str] = None  # 模式1：从FASTQ开始|Mode 1: Start from FASTQ
    bam_file: Optional[str] = None  # 模式2：从BAM开始|Mode 2: Start from BAM (Juicer format)
    hic_file: Optional[str] = None  # 模式3：从.hic开始|Mode 3: Start from .hic file

    # 可选参数|Optional parameters
    threads: int = 64  # 适合300GB内存|Suitable for 300GB memory

    # Juicer参数|Juicer parameters
    restriction_enzyme: str = "MboI"  # 限制性内切酶|Restriction enzyme
    quality_cutoff: int = 30  # MAPQ质量阈值|MAPQ quality cutoff (Juicer default: 30)

    # PlotHiC参数|PlotHiC parameters
    resolution: int = 100000  # 热图分辨率|Heatmap resolution (100kb default for whole genome)
    color_map: str = "YlOrRd"  # 颜色方案|Color scheme (PlotHiC default)
    dpi: int = 300  # 图像分辨率|Image resolution
    output_format: str = "pdf"  # 输出格式|Output format (pdf, png, svg, etc.)

    # 流程控制参数|Process control parameters
    skip_existing: bool = True
    force_juicer: bool = False  # 强制重新运行Juicer|Force rerun Juicer

    # 工具路径|Tool paths
    juicer_dir: str = field(
        default_factory=lambda: get_tool_path(
            'juicer',
            '~/software/juicer',
            'JUICER_DIR'
        )
    )
    bwa_path: str = field(
        default_factory=lambda: get_tool_path(
            'bwa',
            '~/.local/bin/bwa',
            'BWA_PATH'
        )
    )
    samtools_path: str = field(
        default_factory=lambda: get_tool_path(
            'samtools',
            '~/.local/bin/samtools',
            'SAMTOOLS_PATH'
        )
    )
    java_path: str = field(
        default_factory=lambda: get_tool_path(
            'java',
            '~/miniforge3/envs/juicer_v.1.6/bin/java',
            'JAVA_BIN_PATH'
        )
    )
    plothic_path: str = field(
        default_factory=lambda: get_tool_path(
            'plothic',
            '~/miniforge3/envs/plothic_v.1.0.0/bin/plothic',
            'PLOTHIC_PATH'
        )
    )

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开工具路径|Expand tool paths
        self.juicer_dir = Path(expand_path(self.juicer_dir))
        self.bwa_path = expand_path(self.bwa_path)
        self.samtools_path = expand_path(self.samtools_path)
        self.java_path = expand_path(self.java_path)
        self.plothic_path = expand_path(self.plothic_path)

        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 创建Juicer输出目录|Create Juicer output directories
        self.juicer_aligned_dir = self.output_path / "aligned"
        self.juicer_aligned_dir.mkdir(exist_ok=True)
        self.plot_dir = self.output_path / "plot"
        self.plot_dir.mkdir(exist_ok=True)

        # Juicer脚本路径|Juicer script paths
        self.juicer_sh = self.juicer_dir / "scripts" / "juicer.sh"

        # 展开输入文件路径|Expand input file paths
        if self.bam_file:
            self.bam_file = expand_path(self.bam_file)
        if self.hic_file:
            self.hic_file = expand_path(self.hic_file)
        if self.fastq_r1:
            self.fastq_r1 = expand_path(self.fastq_r1)
        if self.fastq_r2:
            self.fastq_r2 = expand_path(self.fastq_r2)

        # 自动提取genome_id（如果未提供）|Auto-extract genome_id if not provided
        if not self.genome_id:
            self.genome_id = self._extract_genome_id()

    def _extract_genome_id(self) -> str:
        """从基因组文件路径提取genome_id|Extract genome_id from genome file path

        Returns:
            str: 基因组ID|Genome ID
        """
        genome_path = Path(self.genome)

        # 尝试从文件名提取|Try to extract from filename
        # 例如: hg19.fa -> hg19, genome.fa -> genome
        stem = genome_path.stem  # 去掉扩展名|Remove extension
        if stem.endswith('.fa'):
            stem = stem[:-3]

        # 如果文件名是genome或assembly，使用目录名
        # If filename is 'genome' or 'assembly', use directory name
        if stem in ['genome', 'assembly']:
            stem = genome_path.parent.name

        return stem

    def get_mode(self):
        """获取运行模式|Get run mode

        Returns:
            str: 'fastq', 'bam', 或 'hic'|'fastq', 'bam', or 'hic'
        """
        if self.hic_file:
            return 'hic'
        elif self.bam_file:
            return 'bam'
        elif self.fastq_r1 and self.fastq_r2:
            return 'fastq'
        else:
            raise ValueError("必须提供输入文件：fastq_r1+fastq_r2, bam_file, 或 hic_file|Must provide input files: fastq_r1+fastq_r2, bam_file, or hic_file")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查基因组文件|Check genome file
        if not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome}")

        # 根据模式检查输入文件|Check input files based on mode
        mode = self.get_mode()

        if mode == 'fastq':
            if not os.path.exists(self.fastq_r1):
                errors.append(f"R1测序文件不存在|R1 fastq file not found: {self.fastq_r1}")
            if not os.path.exists(self.fastq_r2):
                errors.append(f"R2测序文件不存在|R2 fastq file not found: {self.fastq_r2}")
        elif mode == 'bam':
            if not os.path.exists(self.bam_file):
                errors.append(f"BAM文件不存在|BAM file not found: {self.bam_file}")
        elif mode == 'hic':
            if not os.path.exists(self.hic_file):
                errors.append(f"Hi-C文件不存在|Hi-C file not found: {self.hic_file}")

        # 对于非hic模式，检查Juicer脚本|For non-hic modes, check Juicer scripts
        if mode != 'hic':
            if not os.path.exists(self.juicer_sh):
                errors.append(f"Juicer脚本不存在|Juicer script not found: {self.juicer_sh}")

        # 检查参数合法性|Check parameter validity
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if self.resolution <= 0:
            errors.append(f"分辨率必须为正数|Resolution must be positive: {self.resolution}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class HiCProConfig:
    """HiCPro配置类|HiCPro configuration class (HiCPro + PlotHiC)"""

    # 必需参数|Required parameters
    genome: str  # 基因组FASTA文件|Genome FASTA file
    output_dir: str  # 输出目录|Output directory
    genome_id: str  # 基因组ID|Genome ID

    # 输入文件|Input files
    fastq_r1: str  # R1 FASTQ文件
    fastq_r2: str  # R2 FASTQ文件

    # 可选参数|Optional parameters
    threads: int = 64
    max_memory_gb: int = 200  # 最大内存限制（GB）|Maximum memory limit in GB
    restriction_enzyme: str = "MboI"  # 限制性内切酶|Restriction enzyme
    bowtie2_idx: Optional[str] = None  # BWA索引路径|Bowtie2 index path (auto-built if None)

    # Contact map参数|Contact map parameters
    bin_sizes: str = "20000 40000 150000 500000 1000000"  # bin大小列表|Bin sizes
    matrix_format: str = "upper"  # 矩阵格式|Matrix format (upper/lower/complete)

    # PlotHiC参数|PlotHiC parameters
    resolution: int = 100000  # 热图分辨率|Heatmap resolution (100kb default)
    color_map: str = "YlOrRd"  # 颜色方案|Color scheme (PlotHiC default)
    dpi: int = 300  # 图像分辨率|Image resolution
    output_format: str = "pdf"  # 输出格式|Output format (pdf, png, svg, etc.)
    bar_max: int = 1  # 颜色条最大值|Color bar maximum value

    # 流程控制参数|Process control parameters
    skip_existing: bool = True
    force_hicpro: bool = False

    # 工具路径|Tool paths
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
    singularity_exec: str = field(
        default_factory=lambda: get_tool_path(
            'singularity_exec',
            '~/miniforge3/envs/singularity_v.3.8.7/bin/singularity',
            'SINGULARITY_EXEC'
        )
    )
    plothic_path: str = field(
        default_factory=lambda: get_tool_path(
            'plothic',
            '~/miniforge3/envs/plothic_v.1.0.0/bin/plothic',
            'PLOTHIC_PATH'
        )
    )

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开工具路径|Expand tool paths
        self.hicpro_path = expand_path(self.hicpro_path)
        self.plothic_path = expand_path(self.plothic_path)

        # hicpro_sif可选（如果为空则直接用HiC-Pro）|hicpro_sif is optional (use HiC-Pro directly if empty)
        if self.hicpro_sif:
            self.hicpro_sif = expand_path(self.hicpro_sif)
            self.use_singularity = True
        else:
            self.hicpro_sif = None
            self.use_singularity = False

        self.singularity_exec = expand_path(self.singularity_exec)

        # 展开文件路径|Expand file paths
        self.genome = expand_path(self.genome)
        self.fastq_r1 = expand_path(self.fastq_r1)
        self.fastq_r2 = expand_path(self.fastq_r2)

        # 提取样本名|Extract sample name from fastq filename
        r1_name = Path(self.fastq_r1).stem  # 去掉.gz|Remove .gz
        if r1_name.endswith('.fastq') or r1_name.endswith('.fq'):
            r1_name = Path(r1_name).stem  # 再去掉.fastq或.fq|Remove .fastq or .fq
        self.sample_name = (r1_name
                           .replace('_R1', '').replace('_R2', '')
                           .replace('_1', '').replace('_2', '')
                           .replace('.clean', ''))

        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 创建输出目录|Create output directories
        self.hicpro_output_dir = self.output_path / "hicpro_output"
        self.plot_dir = self.output_path / "plot"
        self.plot_dir.mkdir(exist_ok=True)

        # HiCPro配置文件路径|HiCPro config file path
        self.config_file = self.output_path / "hicpro.conf"

        # 热图输出文件|Heatmap output file
        # 使用genome_id作为前缀|Use genome_id as prefix
        self.heatmap_file = self.plot_dir / f"{self.genome_id}_hic_heatmap.{self.output_format}"

        # 如果未提供bowtie2索引，自动生成路径|Auto-generate bowtie2 index path if not provided
        # HiC-Pro使用: BOWTIE2_IDX=${BOWTIE2_IDX_PATH}/${REFERENCE_GENOME}|HiC-Pro uses: BOWTIE2_IDX=${BOWTIE2_IDX_PATH}/${REFERENCE_GENOME}
        # 所以bowtie2_idx应该是包含基因组ID的完整索引前缀|So bowtie2_idx should be full index prefix with genome_id
        if not self.bowtie2_idx:
            # 尝试查找索引|Try to find index
            genome_path = Path(self.genome)

            # 尝试多种可能的索引位置|Try multiple possible index locations
            possible_indices = [
                # 1. genome目录下以genome_id命名的索引|Index named with genome_id in genome dir
                genome_path.parent / self.genome_id,
                # 2. genome目录下以genome_id.fa命名的索引|Index named genome_id.fa in genome dir
                genome_path.parent / f"{self.genome_id}.fa",
                # 3. genome文件去掉.fa的索引|Index from genome file without .fa
                genome_path.parent / genome_path.stem if genome_path.suffix == '.fa' else genome_path,
                # 4. genome文件本身|genome file itself
                genome_path
            ]

            # 找到第一个存在的索引|Find first existing index
            for idx_path in possible_indices:
                # 检查bowtie2索引文件是否存在（只检查.bt2）|Check if bowtie2 index files exist (.bt2 only)
                has_index = False
                if (Path(str(idx_path) + '.1.bt2')).exists():
                    has_index = True
                    break

            if has_index:
                self.bowtie2_idx = str(idx_path)
            else:
                # 没有找到索引，使用genome_id（用户需要手动构建）|No index found, use genome_id (user needs to build manually)
                self.bowtie2_idx = str(genome_path.parent / self.genome_id)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input files
        if not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome}")
        if not os.path.exists(self.fastq_r1):
            errors.append(f"R1测序文件不存在|R1 fastq not found: {self.fastq_r1}")
        if not os.path.exists(self.fastq_r2):
            errors.append(f"R2测序文件不存在|R2 fastq not found: {self.fastq_r2}")

        # 检查HiC-Pro|Check HiC-Pro
        if not os.path.exists(self.hicpro_path):
            errors.append(f"HiC-Pro不存在|HiC-Pro not found: {self.hicpro_path}")

        # 检查singularity镜像（仅在use_singularity时）|Check singularity image (only when use_singularity)
        if self.use_singularity and not os.path.exists(self.hicpro_sif):
            errors.append(f"HiCPro镜像不存在|HiCPro image not found: {self.hicpro_sif}")

        # 检查plothic|Check plothic
        if not os.path.exists(self.plothic_path):
            errors.append(f"plothic不存在|plothic not found: {self.plothic_path}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

