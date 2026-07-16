"""
K-mer工具配置管理模块|K-mer Tools Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path, get_tool_path


@dataclass
class KmerToolsConfig:
    """K-mer工具基础配置类|K-mer Tools Base Configuration Class"""

    # 通用参数|Common parameters
    threads: int = 64  # 线程数|Thread count
    verbose: bool = False  # 详细输出|Verbose output

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        pass


@dataclass
class BuildConfig(KmerToolsConfig):
    """构建k-mer库配置类|Build K-mer Database Configuration Class"""

    # 必需参数|Required parameters
    input_dir: str = ""  # 输入目录（FASTQ/FASTA文件）|Input directory (FASTQ/FASTA files)
    output_dir: str = ""  # 输出目录|Output directory

    # 模式选择|Mode selection
    use_kmtricks: bool = True  # 使用kmtricks模式|Use kmtricks mode (default: True)

    # k-mer参数|K-mer parameters
    kmer_size: int = 51  # k-mer大小|K-mer size
    hard_min: int = 2  # 最小丰度|Minimum abundance
    recurrence_min: int = 1  # 最小重现次数|Minimum recurrence
    header_db_key: str = "kmer_header"  # 数据库中的header key|Header key in database

    # kmtricks管道参数|kmtricks pipeline parameters
    mode: str = "kmer:pa:bin"  # kmtricks模式|kmtricks mode
    minimizer_size: int = 10  # minimizer大小|Minimizer size
    batch_size: int = 20000  # 批量写入大小|Batch write size (for RocksDB import)
    bloom_bits: int = 15  # Bloom filter位数|Bloom filter bits per key
    tmp_dir: str = ""  # kmtricks临时目录|kmtricks temporary directory (default: output_dir/tmp)
    nb_partitions: int = 0  # 分区数，0表示自动计算，-1表示使用kmtricks默认|Partitions, 0=auto, -1=kmtricks default

    # kmindex参数|kmindex parameters
    index_name: str = ""  # 索引名称|Index name (for kmindex)
    bloom_size: int = 1000000000000  # 布隆过滤器大小|Bloom filter size (for kmindex)

    # 工具路径|Tool paths (None时由__post_init__经get_tool_path按优先级解析: 环境变量>配置文件>默认|None -> resolved by get_tool_path in __post_init__: env>config>default)
    kmtricks_path: Optional[str] = None
    kmindex_path: Optional[str] = None
    bgzip_path: Optional[str] = None

    # FOF参数|FOF parameters
    fof_suffix_1: str = "_1.clean.fq.gz"  # R1后缀|R1 suffix
    fof_suffix_2: str = "_2.clean.fq.gz"  # R2后缀|R2 suffix
    fof_file: str = ""  # 预存的FOF文件路径（如果已存在）|Pre-existing FOF file path (if exists)
    header_file: str = ""  # 预存的header文件路径|Pre-existing header file path

    # 内部使用|Internal use
    run_dir: str = ""  # 运行目录|Run directory
    _num_samples: int = 0  # 样本数量|Number of samples (internal)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        super().__post_init__()
        # 工具路径: None时按优先级解析(环境变量>配置文件>默认),再展开~|Tool paths: resolve if None (env>config>default), then expand ~
        if not self.kmtricks_path:
            self.kmtricks_path = get_tool_path('kmtricks', '~/miniforge3/envs/biopytools/bin/kmtricks', 'KMTRICKS_PATH')
        if not self.kmindex_path:
            self.kmindex_path = get_tool_path('kmindex', '~/miniforge3/envs/kmindex_v.0.6.0/bin/kmindex', 'KMINDEX_PATH')
        if not self.bgzip_path:
            self.bgzip_path = get_tool_path('bgzip', 'bgzip', 'BGZIP_PATH')
        self.kmtricks_path = expand_path(self.kmtricks_path)
        self.kmindex_path = expand_path(self.kmindex_path)
        self.bgzip_path = expand_path(self.bgzip_path)
        # 展开用户路径中的~和环境变量|Expand ~ and env vars in user paths
        if self.input_dir:
            self.input_dir = expand_path(self.input_dir)
            self.input_path = Path(self.input_dir)
        if self.output_dir:
            self.output_dir = expand_path(self.output_dir)
            self.output_path = Path(self.output_dir)
            # 创建输出目录|Create output directory
            self.output_path.mkdir(parents=True, exist_ok=True)
            self.output_dir = str(self.output_path.absolute())
            # 设置RocksDB路径|Set RocksDB path
            self.rocksdb_path = self.output_path / "rocksdb"
        # 展开可选输入路径|Expand optional input paths
        if self.fof_file:
            self.fof_file = expand_path(self.fof_file)
        if self.header_file:
            self.header_file = expand_path(self.header_file)
        # 设置默认临时目录|Set default temporary directory
        if not self.tmp_dir and self.output_dir:
            self.tmp_dir = str(Path(self.output_dir) / "tmp")
        if self.tmp_dir:
            self.tmp_dir = expand_path(self.tmp_dir)
            # 创建临时目录|Create temporary directory
            Path(self.tmp_dir).mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入目录|Check input directory
        if not self.input_dir:
            errors.append("输入目录不能为空|Input directory cannot be empty")
        elif not Path(self.input_dir).exists():
            errors.append(f"输入目录不存在|Input directory does not exist: {self.input_dir}")

        # 检查输出目录|Check output directory
        if not self.output_dir:
            errors.append("输出目录不能为空|Output directory cannot be empty")

        # 检查k-mer大小|Check k-mer size
        if self.kmer_size < 8 or self.kmer_size > 127:
            errors.append(f"k-mer大小必须在8-127之间|K-mer size must be between 8-127: {self.kmer_size}")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def calculate_nb_partitions(self, num_samples: int) -> int:
        """根据样本数动态计算最佳分区数|Calculate optimal number of partitions based on sample count

        Args:
            num_samples: 样本数量|Number of samples

        Returns:
            int: 最佳分区数|Optimal number of partitions
        """
        import math

        if num_samples <= 0:
            return 0

        # 基础策略：分区数应该接近样本数，但不要超过线程数太多
        # 公式：min(样本数, 线程数 * 1.5) 向上取整到最近的10的倍数
        base_partitions = min(num_samples, int(self.threads * 1.5))

        # 向上取整到最近的10的倍数
        optimal = math.ceil(base_partitions / 10) * 10

        # 确保至少为10，最多为1000（kmtricks限制）
        optimal = max(10, min(optimal, 1000))

        return optimal


@dataclass
class QueryConfig(KmerToolsConfig):
    """查询k-mer库配置类|Query K-mer Database Configuration Class"""

    # 必需参数|Required parameters
    rocksdb_dir: str = ""  # RocksDB数据库目录|RocksDB database directory
    query_fasta: str = ""  # 查询FASTA文件|Query FASTA file
    output_dir: str = ""  # 输出目录|Output directory
    index_dir: str = ""  # kmindex索引目录|kmindex index directory

    # 模式选择|Mode selection
    use_kmtricks: bool = True  # 使用kmtricks模式|Use kmtricks mode (default: True)

    # 查询参数|Query parameters
    kmer_size: int = 51  # k-mer大小（必须与建库时一致）|K-mer size (must match build)
    header_db_key: str = "kmer_header"  # 数据库中的header key|Header key in database
    bed_file: str = ""  # BED文件路径（用于生成位置丰度文件）|BED file path (for generating position-abundance file)

    # kmindex查询参数|kmindex query parameters
    index_name: str = ""  # 索引名称|Index name (for kmindex)
    zvalue: int = 0  # findere算法z值|findere z-value
    threshold: float = 0.0  # 共享k-mer阈值|Shared k-mer threshold
    output_format: str = "matrix"  # 输出格式|Output format (json|matrix)

    # 工具路径|Tool paths
    bam2fastq_path: str = "bam2fastq"  # 这里实际不会用到，保留用于一致性|Not used, kept for consistency
    kmindex_path: Optional[str] = None  # kmindex路径|kmindex path (None->get_tool_path解析|resolved via get_tool_path)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        super().__post_init__()
        # 工具路径: None时按优先级解析(环境变量>配置文件>默认),再展开~|Tool path: resolve if None (env>config>default), then expand ~
        if not self.kmindex_path:
            self.kmindex_path = get_tool_path('kmindex', '~/miniforge3/envs/kmindex_v.0.6.0/bin/kmindex', 'KMINDEX_PATH')
        self.kmindex_path = expand_path(self.kmindex_path)
        if self.rocksdb_dir:
            self.rocksdb_dir = expand_path(self.rocksdb_dir)
            self.rocksdb_path = Path(self.rocksdb_dir)
        if self.query_fasta:
            self.query_fasta = expand_path(self.query_fasta)
            self.query_fasta_path = Path(self.query_fasta)
        if self.index_dir:
            self.index_dir = expand_path(self.index_dir)
            self.index_path = Path(self.index_dir)
        if self.bed_file:
            self.bed_file = expand_path(self.bed_file)
            self.bed_file_path = Path(self.bed_file)
        if self.output_dir:
            self.output_dir = expand_path(self.output_dir)
            self.output_path = Path(self.output_dir)
            # 创建输出目录|Create output directory
            self.output_path.mkdir(parents=True, exist_ok=True)
            self.output_dir = str(self.output_path.absolute())

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if self.use_kmtricks:
            # 检查RocksDB目录|Check RocksDB directory
            if not self.rocksdb_dir:
                errors.append("RocksDB目录不能为空|RocksDB directory cannot be empty")
            elif not Path(self.rocksdb_dir).exists():
                errors.append(f"RocksDB目录不存在|RocksDB directory does not exist: {self.rocksdb_dir}")
        else:
            # 检查kmindex索引目录|Check kmindex index directory
            if not self.index_dir:
                errors.append("kmindex索引目录不能为空|kmindex index directory cannot be empty")
            elif not Path(self.index_dir).exists():
                errors.append(f"kmindex索引目录不存在|kmindex index directory does not exist: {self.index_dir}")

        # 检查查询FASTA|Check query FASTA
        if not self.query_fasta:
            errors.append("查询FASTA文件不能为空|Query FASTA file cannot be empty")
        elif not Path(self.query_fasta).exists():
            errors.append(f"查询FASTA文件不存在|Query FASTA file does not exist: {self.query_fasta}")

        # 检查输出目录|Check output directory
        if not self.output_dir:
            errors.append("输出目录不能为空|Output directory cannot be empty")

        # 检查k-mer大小|Check k-mer size
        if self.kmer_size <= 0:
            errors.append(f"k-mer大小必须为正数|K-mer size must be positive: {self.kmer_size}")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class SplitFastaConfig(KmerToolsConfig):
    """分割FASTA配置类|Split FASTA Configuration Class"""

    # 必需参数|Required parameters
    input_fasta: str = ""  # 输入FASTA文件|Input FASTA file
    output_dir: str = ""  # 输出目录|Output directory

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        super().__post_init__()
        if self.input_fasta:
            self.input_fasta = expand_path(self.input_fasta)
            self.input_fasta_path = Path(self.input_fasta)
        if self.output_dir:
            self.output_dir = expand_path(self.output_dir)
            self.output_path = Path(self.output_dir)
            # 创建输出目录|Create output directory
            self.output_path.mkdir(parents=True, exist_ok=True)
            self.output_dir = str(self.output_path.absolute())

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入FASTA|Check input FASTA
        if not self.input_fasta:
            errors.append("输入FASTA文件不能为空|Input FASTA file cannot be empty")
        elif not Path(self.input_fasta).exists():
            errors.append(f"输入FASTA文件不存在|Input FASTA file does not exist: {self.input_fasta}")

        # 检查输出目录|Check output directory
        if not self.output_dir:
            errors.append("输出目录不能为空|Output directory cannot be empty")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class GenFofConfig(KmerToolsConfig):
    """生成FOF配置类|Generate FOF Configuration Class"""

    # 必需参数|Required parameters
    input_dir: str = ""  # 输入目录|Input directory
    output_file: str = ""  # 输出FOF文件|Output FOF file

    # FOF参数|FOF parameters
    suffix_1: str = "_1.clean.fq.gz"  # R1后缀|R1 suffix
    suffix_2: str = "_2.clean.fq.gz"  # R2后缀|R2 suffix

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        super().__post_init__()
        if self.input_dir:
            self.input_dir = expand_path(self.input_dir)
            self.input_path = Path(self.input_dir)
        if self.output_file:
            self.output_file = expand_path(self.output_file)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入目录|Check input directory
        if not self.input_dir:
            errors.append("输入目录不能为空|Input directory cannot be empty")
        elif not Path(self.input_dir).exists():
            errors.append(f"输入目录不存在|Input directory does not exist: {self.input_dir}")

        # 检查输出文件|Check output file
        if not self.output_file:
            errors.append("输出FOF文件不能为空|Output FOF file cannot be empty")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class ImportDBConfig(KmerToolsConfig):
    """导入数据库配置类|Import Database Configuration Class"""

    # 必需参数|Required parameters
    input_matrix: str = ""  # 输入矩阵文件（可以是gzip压缩）|Input matrix file (can be gzipped)
    output_db: str = ""  # 输出RocksDB目录|Output RocksDB directory

    # 导入参数|Import parameters
    input_delimiter: str = "\t"  # 输入分隔符|Input delimiter
    batch_size: int = 20000  # 批量写入大小|Batch write size
    bloom_bits: int = 15  # Bloom filter位数|Bloom filter bits per key
    force_overwrite: bool = False  # 强制覆盖|Force overwrite

    # Header参数|Header parameters
    header_file: Optional[str] = None  # Header文件路径|Header file path
    header_db_key: str = "kmer_header"  # 数据库中的header key|Header key in database

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        super().__post_init__()
        if self.input_matrix:
            self.input_matrix = expand_path(self.input_matrix)
            self.input_matrix_path = Path(self.input_matrix)
        if self.output_db:
            self.output_db = expand_path(self.output_db)
            self.output_db_path = Path(self.output_db)
            # 创建输出目录|Create output directory
            self.output_db_path.parent.mkdir(parents=True, exist_ok=True)
            self.output_db = str(self.output_db_path.absolute())
        if self.header_file:
            self.header_file = expand_path(self.header_file)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入矩阵|Check input matrix
        if not self.input_matrix:
            errors.append("输入矩阵文件不能为空|Input matrix file cannot be empty")
        elif not Path(self.input_matrix).exists():
            errors.append(f"输入矩阵文件不存在|Input matrix file does not exist: {self.input_matrix}")

        # 检查输出数据库|Check output database
        if not self.output_db:
            errors.append("输出RocksDB目录不能为空|Output RocksDB directory cannot be empty")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        # 检查批量大小|Check batch size
        if self.batch_size <= 0:
            errors.append(f"批量写入大小必须为正数|Batch write size must be positive: {self.batch_size}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


@dataclass
class ExtractConfig(KmerToolsConfig):
    """提取k-mer配置类|Extract K-mer Configuration Class"""

    # 必需参数|Required parameters
    fasta_file: str = ""  # 输入FASTA文件|Input FASTA file
    output_dir: str = ""  # 输出目录|Output directory

    # k-mer参数|K-mer parameters
    kmer_size: int = 51  # k-mer大小|K-mer size

    # 提取方法|Extraction method
    extract_method: str = "unikmer"  # 提取方法：unikmer或pyfastx|Extraction method: unikmer or pyfastx
    unikmer_path: Optional[str] = None  # unikmer可执行文件路径|unikmer executable path (None->get_tool_path解析|resolved via get_tool_path)

    # 输出文件|Output files
    kmer_output: str = ""  # kmer FASTA文件|Kmer FASTA file
    kmer_pos_output: str = ""  # kmer位置文件|Kmer position file
    extract_output_bed: bool = True  # 是否输出BED文件|Whether to output BED file

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        super().__post_init__()
        # 工具路径: None时按优先级解析(环境变量>配置文件>默认),再展开~|Tool path: resolve if None (env>config>default), then expand ~
        if not self.unikmer_path:
            self.unikmer_path = get_tool_path('unikmer', 'unikmer', 'UNIKMER_PATH')
        self.unikmer_path = expand_path(self.unikmer_path)
        if self.fasta_file:
            self.fasta_file = expand_path(self.fasta_file)
            self.fasta_path = Path(self.fasta_file)
        if self.output_dir:
            self.output_dir = expand_path(self.output_dir)
            self.output_path = Path(self.output_dir)
            # 创建输出目录|Create output directory
            self.output_path.mkdir(parents=True, exist_ok=True)
            self.output_dir = str(self.output_path.absolute())

        # 如果没有指定kmer_output，使用默认值|If kmer_output not specified, use default
        if not self.kmer_output and self.fasta_file and self.output_dir:
            base_name = Path(self.fasta_file).stem
            self.kmer_output = str(Path(self.output_dir) / f"{base_name}_kmer_{self.kmer_size}.fa")

        # 如果没有指定kmer_pos_output，使用默认值|If kmer_pos_output not specified, use default
        if not self.kmer_pos_output and self.fasta_file and self.output_dir:
            base_name = Path(self.fasta_file).stem
            self.kmer_pos_output = str(Path(self.output_dir) / f"{base_name}_kmer_{self.kmer_size}_pos.txt")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入FASTA|Check input FASTA
        if not self.fasta_file:
            errors.append("输入FASTA文件不能为空|Input FASTA file cannot be empty")
        elif not Path(self.fasta_file).exists():
            errors.append(f"输入FASTA文件不存在|Input FASTA file does not exist: {self.fasta_file}")

        # 检查输出目录|Check output directory
        if not self.output_dir:
            errors.append("输出目录不能为空|Output directory cannot be empty")

        # 检查k-mer大小|Check k-mer size
        if self.kmer_size <= 0:
            errors.append(f"k-mer大小必须为正数|K-mer size must be positive: {self.kmer_size}")

        # 检查提取方法|Check extraction method
        if self.extract_method not in ["unikmer", "pyfastx"]:
            errors.append(f"提取方法必须是unikmer或pyfastx|Extraction method must be unikmer or pyfastx: {self.extract_method}")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
