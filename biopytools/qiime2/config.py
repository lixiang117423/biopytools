"""
QIIME2配置管理模块|QIIME2 Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path, get_tool_path


# 扩增子类型|Amplicon types
VALID_AMPLICON = {'16s', 'its'}

# 聚类方法|Clustering methods
VALID_METHOD = {'asv', 'otu'}

# 默认V3-V4引物(341F/806R)|Default V3-V4 primers (341F/806R)
DEFAULT_FWD_PRIMER = 'CCTACGGGNGGCWGCAG'   # 341F
DEFAULT_REV_PRIMER = 'GACTACHVGGGTATCTAATCC'  # 806R

# IUPAC核苷酸编码|IUPAC nucleotide codes
IUPAC_CODES = set('ACGTURYSWKMBDHVNacgturyswkmbdhvn')

# 默认数据库与工具路径|Default database and tool paths
DEFAULT_DATABASE_DIR = '~/database/qiime2'
DEFAULT_QIIME_BIN = '~/miniforge3/envs/qiime_v.2024.10.1/bin/qiime'

# SILVA参考文件名|SILVA reference filenames
SILVA_NR99_FNAME = 'SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz'

# UNITE参考文件名(99%阈值)|UNITE reference filenames (99% threshold)
UNITE_TGZ_FNAME = 'UNITE_qiime_release_19.02.2025.tgz'
UNITE_FASTA_FNAME = 'sh_refs_qiime_ver10_99_19.02.2025.fasta'
UNITE_TAXONOMY_FNAME = 'sh_taxonomy_qiime_ver10_99_19.02.2025.txt'


@dataclass
class Qiime2Config:
    """QIIME2分析配置类|QIIME2 Analysis Configuration Class"""

    # ===== 必需参数|Required parameters =====
    input_dir: str

    # ===== 输出配置|Output configuration =====
    output_dir: str = './qiime2_output'

    # ===== 分析类型|Analysis type =====
    amplicon: str = '16s'    # 16s 或 its|16s or its
    method: str = 'asv'      # asv (DADA2) 或 otu (vsearch)|asv (DADA2) or otu (vsearch)

    # ===== 引物与截断|Primers and truncation =====
    fwd_primer: str = DEFAULT_FWD_PRIMER
    rev_primer: str = DEFAULT_REV_PRIMER
    trunc_len_f: int = 0     # 0表示不截断|0 means no truncation
    trunc_len_r: int = 0     # 0表示不截断|0 means no truncation
    trim_left_f: int = 0
    trim_left_r: int = 0

    # ===== 抽样深度|Sampling depth =====
    # 0=自动(取每样本reads总数的第10百分位)|0=auto (10th percentile of per-sample totals)
    sampling_depth: int = 0

    # ===== OTU聚类/分类参数|OTU clustering / classification params =====
    perc_identity: float = 0.97   # OTU聚类相似度|OTU clustering identity
    confidence: float = 0.7       # classify-sklearn置信度|classify-sklearn confidence
    min_length: int = 50          # extract-reads最小长度|extract-reads min length
    max_length: int = 0           # 0=不限制|0=no limit

    # ===== 运行控制|Run control =====
    threads: int = 12
    validate_level: str = 'min'   # tools import校验级别|import validate level (min更快|faster)

    # ===== 分类器配置|Classifier configuration =====
    classifier: Optional[str] = None  # 预训练分类器路径|Pre-trained classifier path (None=自动训练|auto-train)

    # ===== 路径配置|Path configuration (支持~展开|supports ~) =====
    database_dir: str = field(
        default_factory=lambda: get_tool_path('qiime2_db', DEFAULT_DATABASE_DIR, 'QIIME2_DB')
    )
    qiime_bin: str = field(
        default_factory=lambda: get_tool_path('qiime', DEFAULT_QIIME_BIN, 'QIIME_PATH')
    )
    classifier_cache_dir: Optional[str] = None  # 默认<database_dir>/classifier_cache|default

    # ===== 样品命名|Sample naming =====
    r1_suffix: str = '_1.clean.fq.gz'
    r2_suffix: str = '_2.clean.fq.gz'

    # ===== 跳过控制|Skip control =====
    skip_cutadapt: bool = False   # 跳过引物切除(数据已去引物)|skip primer trimming
    skip_phylogeny: bool = False  # 跳过系统发育(ITS自动True)|skip phylogeny (auto True for ITS)

    # ===== 元数据与杂项|Metadata and misc =====
    metadata_file: Optional[str] = None
    force: bool = False
    verbose: bool = False

    # ===== 内部路径(__post_init__填充)|Internal paths (populated in __post_init__) =====
    output_path: Optional[Path] = field(default=None, init=False, repr=False)
    info_dir: Optional[str] = field(default=None, init=False, repr=False)
    import_dir: Optional[str] = field(default=None, init=False, repr=False)
    trim_dir: Optional[str] = field(default=None, init=False, repr=False)
    denoise_dir: Optional[str] = field(default=None, init=False, repr=False)
    otu_dir: Optional[str] = field(default=None, init=False, repr=False)
    taxonomy_dir: Optional[str] = field(default=None, init=False, repr=False)
    phylogeny_dir: Optional[str] = field(default=None, init=False, repr=False)
    diversity_dir: Optional[str] = field(default=None, init=False, repr=False)
    export_dir: Optional[str] = field(default=None, init=False, repr=False)
    log_dir: Optional[str] = field(default=None, init=False, repr=False)
    work_dir: Optional[str] = field(default=None, init=False, repr=False)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有含~的路径|Expand all paths with ~ (§11.3.1 关键|critical)
        self.input_dir = expand_path(self.input_dir)
        self.output_dir = expand_path(self.output_dir)
        self.database_dir = expand_path(self.database_dir)
        self.qiime_bin = expand_path(self.qiime_bin)
        if self.classifier:
            self.classifier = expand_path(self.classifier)
        if self.metadata_file:
            self.metadata_file = expand_path(self.metadata_file)
        if self.classifier_cache_dir:
            self.classifier_cache_dir = expand_path(self.classifier_cache_dir)
        else:
            # 默认缓存目录置于数据库目录下|Default cache dir under database dir
            self.classifier_cache_dir = os.path.join(self.database_dir, 'classifier_cache')

        # ITS: 系统发育指标无意义,自动跳过|ITS: phylogenetic metrics meaningless, auto-skip
        if self.amplicon == 'its':
            self.skip_phylogeny = True

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准子目录(§12)|Standard subdirectories
        self.info_dir = os.path.join(self.output_dir, '00_pipeline_info')
        self.import_dir = os.path.join(self.output_dir, '01_import')
        self.trim_dir = os.path.join(self.output_dir, '02_trim')
        self.denoise_dir = os.path.join(self.output_dir, '03_denoise')
        self.otu_dir = os.path.join(self.output_dir, '04_otu')
        self.taxonomy_dir = os.path.join(self.output_dir, '05_taxonomy')
        self.phylogeny_dir = os.path.join(self.output_dir, '06_phylogeny')
        self.diversity_dir = os.path.join(self.output_dir, '07_diversity')
        self.export_dir = os.path.join(self.output_dir, '08_export')
        self.log_dir = os.path.join(self.output_dir, '99_logs')
        self.work_dir = os.path.join(self.output_dir, 'work')

        for dir_path in [self.info_dir, self.import_dir, self.trim_dir,
                         self.denoise_dir, self.otu_dir, self.taxonomy_dir,
                         self.phylogeny_dir, self.diversity_dir, self.export_dir,
                         self.log_dir, self.work_dir]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 输入目录|Input directory
        if not os.path.isdir(self.input_dir):
            errors.append(f"输入目录不存在|Input directory not found: {self.input_dir}")

        # 扩增子与方法|Amplicon and method
        if self.amplicon not in VALID_AMPLICON:
            errors.append(f"无效的扩增子类型|Invalid amplicon: {self.amplicon} (可选|options: {VALID_AMPLICON})")

        if self.method not in VALID_METHOD:
            errors.append(f"无效的方法|Invalid method: {self.method} (可选|options: {VALID_METHOD})")

        # 线程数|Threads
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive: {self.threads}")

        # 引物IUPAC校验|Primer IUPAC validation
        for name, primer in [('fwd_primer', self.fwd_primer), ('rev_primer', self.rev_primer)]:
            if not primer:
                errors.append(f"{name}不能为空|{name} cannot be empty")
            else:
                invalid = set(primer) - IUPAC_CODES
                if invalid:
                    errors.append(
                        f"{name}含非法字符|{name} contains invalid characters: {invalid} "
                        f"(仅允许IUPAC编码|only IUPAC codes allowed)"
                    )

        # 截断长度|Truncation lengths
        if self.trunc_len_f < 0 or self.trunc_len_r < 0:
            errors.append(f"截断长度不能为负|Truncation lengths cannot be negative: "
                          f"trunc_len_f={self.trunc_len_f}, trunc_len_r={self.trunc_len_r}")

        # OTU相似度|OTU identity
        if not (0.0 < self.perc_identity <= 1.0):
            errors.append(f"perc_identity必须在(0,1]区间|perc_identity must be in (0,1]: {self.perc_identity}")

        # 置信度|Confidence
        if not (0.0 <= self.confidence <= 1.0):
            errors.append(f"confidence必须在[0,1]区间|confidence must be in [0,1]: {self.confidence}")

        # 抽样深度|Sampling depth
        if self.sampling_depth < 0:
            errors.append(f"sampling_depth不能为负|sampling_depth cannot be negative: {self.sampling_depth}")

        # 预训练分类器|Pre-trained classifier
        if self.classifier and not os.path.exists(self.classifier):
            errors.append(f"分类器文件不存在|Classifier file not found: {self.classifier}")

        # 元数据|Metadata
        if self.metadata_file and not os.path.exists(self.metadata_file):
            errors.append(f"元数据文件不存在|Metadata file not found: {self.metadata_file}")

        # 校验级别|Validate level
        if self.validate_level not in ('min', 'max'):
            errors.append(f"validate_level必须为min或max|validate_level must be min or max: {self.validate_level}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    @property
    def is_16s(self) -> bool:
        """是否为16S分析|Whether 16S analysis"""
        return self.amplicon == '16s'

    @property
    def ref_database_path(self) -> str:
        """当前扩增子对应的原始参考数据库路径|Raw reference DB path for current amplicon"""
        if self.is_16s:
            return os.path.join(self.database_dir, SILVA_NR99_FNAME)
        else:
            return os.path.join(self.database_dir, UNITE_TGZ_FNAME)
