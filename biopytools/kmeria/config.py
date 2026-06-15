"""KMERIA配置管理模块|KMERIA Configuration Management Module"""

import os
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List
from ..common.paths import expand_path


# GWAS场景下k-mer大小上限（GitHub README明确规定）
# Maximum k-mer size for GWAS scenarios (clearly stated in GitHub README)
KMER_GWAS_MAX = 31
# 软件支持的最大k-mer大小（kcount.c校验）
# Maximum k-mer size supported by the software (kcount.c validation)
KMER_SOFT_MAX = 63


@dataclass
class KMERIAConfig:
    """KMERIA基础配置类|KMERIA Base Configuration Class"""

    # 软件路径配置|Software path configuration
    kmeria_path: str = '~/software/kmeria'
    kmeria_conda: str = '~/miniforge3/envs/kmeriaenv'

    # 行为控制|Behavior control
    force: bool = False  # 强制重新运行|Force re-run even if output exists

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.kmeria_path = os.path.normpath(os.path.abspath(expand_path(self.kmeria_path)))
        self.kmeria_conda = os.path.normpath(os.path.abspath(expand_path(self.kmeria_conda)))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not os.path.exists(self.kmeria_path):
            errors.append(f"KMERIA路径不存在|KMERIA path does not exist: {self.kmeria_path}")

        if not os.path.exists(self.kmeria_conda):
            errors.append(f"KMERIA conda环境不存在|KMERIA conda env does not exist: {self.kmeria_conda}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class CountConfig(KMERIAConfig):
    """k-mer计数配置类|K-mer Count Configuration Class"""

    # 必需参数|Required parameters
    fastq_dir: str = ''
    samples_file: str = ''  # 样本列表文件|Sample list file
    output_dir: str = './01_kmer_counts'

    # k-mer参数|K-mer parameters
    kmer_size: int = 31  # KMERIA v2.0.4支持2-63（GWAS建议≤31）|Supports 2-63 (≤31 recommended for GWAS)
    count_separate_strands: bool = False  # -C: count strands separately (no canonical)
    text_output: bool = False  # -T: text output instead of binary

    # 丰度过滤参数|Abundance filter parameters (v2.0.2+)
    min_abund: int = 1  # -m: 最小k-mer丰度|Minimum k-mer abundance (default: 1)
    max_abund: int = -1  # -M: 最大k-mer丰度（-1表示不限）|Maximum k-mer abundance (-1 means unlimited)
    hist_output: Optional[str] = None  # -H: 直方图输出文件|Histogram output file

    # 性能参数|Performance parameters
    threads: int = 24
    batch_size: int = 4  # 每批处理的样本数|Samples per batch

    # 输入文件模式|Input file pattern
    file_pattern: str = '*.fq.gz'  # FASTQ文件模式|FASTQ file pattern

    def __post_init__(self):
        super().__post_init__()
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        if self.fastq_dir:
            self.fastq_dir = os.path.normpath(os.path.abspath(self.fastq_dir))

        if self.samples_file:
            self.samples_file = os.path.normpath(os.path.abspath(self.samples_file))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        super().validate()
        errors = []

        if not self.fastq_dir or not os.path.exists(self.fastq_dir):
            errors.append(f"FASTQ目录不存在|FASTQ directory does not exist: {self.fastq_dir}")

        if not self.samples_file or not os.path.exists(self.samples_file):
            errors.append(f"样本列表文件不存在|Sample list file does not exist: {self.samples_file}")

        if self.kmer_size < 2 or self.kmer_size > KMER_SOFT_MAX:
            errors.append(f"K-mer大小必须在2-{KMER_SOFT_MAX}之间|K-mer size must be between 2-{KMER_SOFT_MAX}: {self.kmer_size}")
        elif self.kmer_size > KMER_GWAS_MAX:
            # 单步count允许>31（用于标准k-mer分析），但下游GWAS不支持
            # Standalone count allows >31 (for standard k-mer analysis), but downstream GWAS is not supported
            logging.getLogger("biopytools.kmeria").warning(
                f"K-mer大小>{KMER_GWAS_MAX}({self.kmer_size})仅适用于标准k-mer分析，下游GWAS流程(m2b/asso)不支持|"
                f"K-mer size >{KMER_GWAS_MAX} ({self.kmer_size}) is for standard k-mer analysis only, "
                f"downstream GWAS pipeline (m2b/asso) is not supported"
            )

        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class KctmConfig(KMERIAConfig):
    """k-mer矩阵构建配置类|K-mer Matrix Construction Configuration Class"""

    # 必需参数|Required parameters
    input_dir: str = ''  # k-mer文件目录|Directory with k-mer files
    output_dir: str = './02_kmer_matrices'

    # 性能参数|Performance parameters
    threads: int = 24
    kmer_batch_size: int = 1000  # -b: 每批k-mer数（与GitHub kmatrix.cpp默认1000一致）|K-mers per batch (matches GitHub kmatrix.cpp default 1000)
    block_size: int = 10000000  # -s: 每块条目数|Entries per block (GitHub default: 10000000)

    def __post_init__(self):
        super().__post_init__()
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        if self.input_dir:
            self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        super().validate()
        errors = []

        if not self.input_dir or not os.path.exists(self.input_dir):
            errors.append(f"输入目录不存在|Input directory does not exist: {self.input_dir}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class FilterConfig(KMERIAConfig):
    """k-mer矩阵过滤配置类|K-mer Matrix Filter Configuration Class"""

    # 必需参数|Required parameters
    input_dir: str = ''  # k-mer矩阵目录|Directory with k-mer matrices
    output_dir: str = './03_filtered_matrices'
    depth_file: str = ''  # 样本深度文件|Sample depth file

    # 过滤参数|Filter parameters (默认值与GitHub kfilter.h一致|Defaults match GitHub kfilter.h)
    max_abund: int = 1000  # -c: 最大k-mer覆盖度|Maximum k-mer coverage (GitHub default: 1000)
    missing_ratio: float = 0.8  # -s: 缺失率阈值|Missing ratio threshold (GitHub default: 0.8, k-mer通常缺失率高)
    ploidy: int = 4  # -p: 基因组倍性|Genome ploidy (GitHub default: 4, 适配多倍体)

    # 性能参数|Performance parameters
    threads: int = 24

    def __post_init__(self):
        super().__post_init__()
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        if self.input_dir:
            self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        if self.depth_file:
            self.depth_file = os.path.normpath(os.path.abspath(self.depth_file))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        super().validate()
        errors = []

        if not self.input_dir or not os.path.exists(self.input_dir):
            errors.append(f"输入目录不存在|Input directory does not exist: {self.input_dir}")

        if not self.depth_file or not os.path.exists(self.depth_file):
            errors.append(f"深度文件不存在|Depth file does not exist: {self.depth_file}")

        if self.missing_ratio < 0 or self.missing_ratio > 1:
            errors.append(f"缺失率必须在0-1之间|Missing ratio must be between 0-1: {self.missing_ratio}")

        if self.ploidy <= 0:
            errors.append(f"倍性必须为正数|Ploidy must be positive: {self.ploidy}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class M2bConfig(KMERIAConfig):
    """k-mer矩阵转BIMBAM格式配置类|Matrix to BIMBAM Configuration Class"""

    # 必需参数|Required parameters
    input_dir: str = ''  # 过滤后的矩阵目录|Directory with filtered matrices
    output_dir: str = './04_bimbam'

    # 归一化参数|Normalization parameters
    normalize: bool = True  # 是否归一化|Whether to normalize
    min_range: float = 0.0  # 最小值|Minimum value
    max_range: float = 2.0  # 最大值|Maximum value
    quantile_norm: bool = False  # 分位数归一化|Quantile normalization
    lower_quantile: float = 0.05  # 下分位数|Lower quantile
    upper_quantile: float = 0.95  # 上分位数|Upper quantile

    # 压缩参数|Compression parameters
    compress: bool = True  # 是否压缩|Whether to compress
    compression_level: int = 6  # 压缩级别|Compression level (0-9)
    bgzf_threads: int = 4  # BGZF线程数|BGZF threads

    # 性能参数|Performance parameters
    threads: int = 24

    def __post_init__(self):
        super().__post_init__()
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        if self.input_dir:
            self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        super().validate()
        errors = []

        if not self.input_dir or not os.path.exists(self.input_dir):
            errors.append(f"输入目录不存在|Input directory does not exist: {self.input_dir}")

        if self.compression_level < 0 or self.compression_level > 9:
            errors.append(f"压缩级别必须在0-9之间|Compression level must be between 0-9: {self.compression_level}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class AssoConfig(KMERIAConfig):
    """k-mer关联分析配置类|K-mer Association Analysis Configuration Class"""

    # 必需参数|Required parameters
    input_dir: str = ''  # BIMBAM文件目录|Directory with BIMBAM files
    pheno_file: str = ''  # 表型文件|Phenotype file
    output_dir: str = './05_association'

    # 分析参数|Analysis parameters
    pheno_col: int = 1  # -n: 表型列|Phenotype column (1-based)
    covar_file: Optional[str] = None  # -c: 协变量文件|Covariate file
    kinship_file: Optional[str] = None  # -k: 亲缘关系矩阵|Kinship matrix
    kinship_precision: int = 10  # --kin-precision: kinship精度|Kinship precision
    output_precision: int = 5  # --out-precision: 输出精度|Output precision

    # 工具选择|Tool selection
    use_bimbam_tools: bool = True  # 使用bimbamAsso|Use bimbamAsso (False=gemma)
    analysis_method: Optional[str] = None  # -m: 分析方法（default/lm/lmm）|Analysis method
    kinship_method: int = 3  # --kin-method: 1=IBS均值, 2=IBS随机, 3=Balding-Nichols|Kinship method (default 3)
    use_kinship: bool = True  # --no-kinship 关闭则不使用kinship|If False, skip kinship
    disable_gls: bool = False  # --disable-gls: 禁用GLS，改用OLS|Disable GLS, use OLS
    write_eigen: bool = False  # --write-eigen: 输出特征值/特征向量|Write eigenvalues/eigenvectors

    # 质量控制|Quality control
    minor_allele_freq: Optional[float] = None  # --maf: 次等位基因频率|Minor allele frequency
    missing_threshold: Optional[float] = None  # --miss: 缺失阈值|Missing threshold

    # 分块分析（大样本）|Block analysis (large samples)
    start_marker: Optional[int] = None  # --start-marker: 起始marker索引|Start marker index
    end_marker: Optional[int] = None  # --end-marker: 结束marker索引|End marker index

    # 其他选项|Other options
    generate_plots: bool = False  # --generate-plots: 生成图表|Generate plots
    compress_output: bool = False  # --compress: 压缩输出|Compress output
    verbose: bool = False  # --verbose: 详细输出|Verbose output
    dry_run: bool = False  # --dry-run: 仅显示命令不执行|Show commands without executing
    no_validate: bool = False  # --no-validate: 跳过输入验证|Skip input validation
    no_check_deps: bool = False  # --no-check-deps: 跳过依赖检查|Skip dependency check
    no_cleanup: bool = False  # --no-cleanup: 保留临时文件|Keep temporary files

    # 性能参数|Performance parameters
    threads: int = 64  # 关联分析线程数|Association threads

    def __post_init__(self):
        super().__post_init__()
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        if self.input_dir:
            self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        if self.pheno_file:
            self.pheno_file = os.path.normpath(os.path.abspath(self.pheno_file))
        if self.covar_file:
            self.covar_file = os.path.normpath(os.path.abspath(self.covar_file))
        if self.kinship_file:
            self.kinship_file = os.path.normpath(os.path.abspath(self.kinship_file))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        super().validate()
        errors = []

        if not self.input_dir or not os.path.exists(self.input_dir):
            errors.append(f"输入目录不存在|Input directory does not exist: {self.input_dir}")

        if not self.pheno_file or not os.path.exists(self.pheno_file):
            errors.append(f"表型文件不存在|Phenotype file does not exist: {self.pheno_file}")

        if self.covar_file and not os.path.exists(self.covar_file):
            errors.append(f"协变量文件不存在|Covariate file does not exist: {self.covar_file}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class PipelineConfig(KMERIAConfig):
    """KMERIA完整流程配置类|KMERIA Complete Pipeline Configuration Class"""

    # 必需参数|Required parameters
    fastq_dir: str = ''  # FASTQ文件目录|FASTQ files directory
    samples_file: str = ''  # 样本列表文件|Sample list file
    depth_file: str = ''  # 测序深度文件|Sequencing depth file
    pheno_file: str = ''  # 表型文件|Phenotype file
    output_dir: str = './kmeria_results'

    # k-mer参数|K-mer parameters
    kmer_size: int = 31
    max_abund: int = 1000

    # 过滤参数|Filter parameters (与GitHub kfilter.h默认值一致|Match GitHub kfilter.h defaults)
    missing_ratio: float = 0.8  # k-mer通常缺失率高，0.8为GitHub默认值|High k-mer missing rate, 0.8 is GitHub default
    ploidy: int = 4  # 默认多倍体，适配甘蔗/甘薯等|Default polyploid for sugarcane/sweetpotato

    # 流程控制|Pipeline control
    step: Optional[str] = None  # 从指定步骤开始|Start from specified step
    steps: List[str] = field(default_factory=list)  # 要运行的步骤|Steps to run

    # 性能参数|Performance parameters
    threads: int = 24
    batch_size: int = 4

    # 关联分析参数|Association parameters
    pheno_col: int = 1
    kinship_file: Optional[str] = None
    covar_file: Optional[str] = None

    # 可选功能|Optional features
    enable_qc: bool = True  # 启用质控统计|Enable QC statistics

    # Post-GWAS参数|Post-GWAS parameters
    genome_file: Optional[str] = None  # 参考基因组|Reference genome
    gff_file: Optional[str] = None  # GFF注释文件|GFF annotation file (可选|optional)
    sample_ratio: float = 0.1  # 高p值位点抽样比例|Sampling ratio for high p-value loci (default: 0.1 = 10%)
    window_size: int = 200000  # 基因查找窗口大小|Gene search window size (default: 200kb)

    # 比对工具选择|Alignment tool selection
    alignment_tool: str = 'bwa'  # 比对工具: 'bwa' (默认|default) 或 'blast' (可选|optional) | Alignment tool: 'bwa' (default) or 'blast' (optional)

    # BWA参数|BWA parameters
    bwa_k: int = 9  # BWA mem -k 参数 (最小种子长度)|Minimum seed length (default: 9)
    bwa_T: int = 10  # BWA mem -T 参数 (最小输出分数)|Minimum score to output (default: 10)
    as_ratio: float = 0.95  # AS (alignment score) 过滤阈值，保留AS >= 最高AS * as_ratio的所有比对|AS filtering threshold, keep alignments with AS >= max_AS * ratio (default: 0.95)

    def __post_init__(self):
        super().__post_init__()
        self.output_path = Path(self.output_dir)

        # 创建各步骤输出目录|Create output directories for each step
        self.dirs = {
            'counts': self.output_path / '01_kmer_counts',
            'matrices': self.output_path / '02_kmer_matrices',
            'filtered': self.output_path / '03_filtered_matrices',
            'bimbam': self.output_path / '04_bimbam',
            'association': self.output_path / '05_association',
            'qc': self.output_path / '06_qc_reports',
            'post_gwas': self.output_path / '07_post_gwas'
        }

        for dir_path in self.dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        if self.fastq_dir:
            self.fastq_dir = os.path.normpath(os.path.abspath(self.fastq_dir))
        if self.samples_file:
            self.samples_file = os.path.normpath(os.path.abspath(self.samples_file))
        if self.depth_file:
            self.depth_file = os.path.normpath(os.path.abspath(self.depth_file))
        if self.pheno_file:
            self.pheno_file = os.path.normpath(os.path.abspath(self.pheno_file))
        if self.kinship_file:
            self.kinship_file = os.path.normpath(os.path.abspath(self.kinship_file))
        if self.covar_file:
            self.covar_file = os.path.normpath(os.path.abspath(self.covar_file))
        if self.genome_file:
            self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))
        if self.gff_file:
            self.gff_file = os.path.normpath(os.path.abspath(self.gff_file))

        # 默认运行所有步骤|Default: run all steps
        if not self.steps:
            self.steps = ['count', 'kctm', 'filter', 'm2b', 'asso']

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        super().validate()
        errors = []

        if not self.fastq_dir or not os.path.exists(self.fastq_dir):
            errors.append(f"FASTQ目录不存在|FASTQ directory does not exist: {self.fastq_dir}")

        if not self.samples_file or not os.path.exists(self.samples_file):
            errors.append(f"样本列表文件不存在|Sample list file does not exist: {self.samples_file}")

        if not self.depth_file or not os.path.exists(self.depth_file):
            errors.append(f"深度文件不存在|Depth file does not exist: {self.depth_file}")

        if not self.pheno_file or not os.path.exists(self.pheno_file):
            errors.append(f"表型文件不存在|Phenotype file does not exist: {self.pheno_file}")

        if self.step and self.step not in ['count', 'kctm', 'filter', 'm2b', 'asso']:
            errors.append(f"无效的步骤名称|Invalid step name: {self.step}")

        # Pipeline必然走GWAS流程，k-mer大小>31会导致m2b/asso失败
        # Pipeline always runs GWAS, k-mer size >31 will break m2b/asso
        if self.kmer_size < 2 or self.kmer_size > KMER_SOFT_MAX:
            errors.append(f"K-mer大小必须在2-{KMER_SOFT_MAX}之间|K-mer size must be between 2-{KMER_SOFT_MAX}: {self.kmer_size}")
        elif self.kmer_size > KMER_GWAS_MAX:
            errors.append(
                f"Pipeline模式下K-mer大小不能超过{KMER_GWAS_MAX}（当前: {self.kmer_size}），"
                f"下游GWAS流程不支持|In pipeline mode, k-mer size must not exceed "
                f"{KMER_GWAS_MAX} (current: {self.kmer_size}), downstream GWAS not supported"
            )

        if self.genome_file and not os.path.exists(self.genome_file):
            errors.append(f"基因组文件不存在|Genome file does not exist: {self.genome_file}")

        if self.gff_file and not os.path.exists(self.gff_file):
            errors.append(f"GFF文件不存在|GFF file does not exist: {self.gff_file}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
