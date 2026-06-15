"""
VG配置类|VG Configuration Classes
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List


def expand_path(path: str) -> str:
    """
    展开路径中的~和环境变量|Expand ~ and environment variables in path

    Args:
        path: 输入路径|Input path

    Returns:
        展开后的绝对路径|Expanded absolute path
    """
    return str(Path(path).expanduser().resolve())


@dataclass
class VGBaseConfig:
    """VG基础配置类|VG Base Configuration Class"""

    # Conda环境|Conda environment
    vg_env: str = "vg_v.1.7.0"  # conda环境名称或路径|conda env name or path

    # 日志|Logging
    log_level: str = "INFO"
    verbose: bool = True

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.log_level = self.log_level.upper()

        # 处理conda环境|Handle conda environment
        if self.vg_env:
            if '/' in self.vg_env or '\\' in self.vg_env:
                # 是路径，提取环境名称|It's a path, extract env name
                env_path = Path(self.vg_env)
                self.vg_env = env_path.name


@dataclass
class VGConstructConfig(VGBaseConfig):
    """VG Construct配置类|VG Construct Configuration Class"""

    # 输入文件|Input files
    reference: str = ""  # 参考基因组FASTA文件|Reference FASTA file
    vcf: str = ""  # VCF文件|VCF file
    output: str = ""  # 输出VG文件|Output VG file

    # 构建参数|Construction parameters
    region: Optional[str] = None  # 指定染色体区域|Specify chromosome region
    threads: int = 12  # 线程数|Number of threads
    node_max: int = 32  # 最大节点大小|Maximum node size
    alt_paths: bool = False  # 保存alt等位基因路径|Save alt allele paths
    progress: bool = True  # 显示进度|Show progress

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        super().__post_init__()

        # 展开路径|Expand paths
        if self.reference:
            self.reference = expand_path(self.reference)
        if self.vcf:
            self.vcf = expand_path(self.vcf)
        if self.output:
            self.output = expand_path(self.output)

        # 确保输出目录存在|Ensure output directory exists
        if self.output:
            Path(self.output).parent.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        if not self.reference:
            raise ValueError("参考基因组文件不能为空|Reference file cannot be empty")

        if not Path(self.reference).exists():
            raise ValueError(f"参考基因组文件不存在|Reference file does not exist: {self.reference}")

        if not self.vcf:
            raise ValueError("VCF文件不能为空|VCF file cannot be empty")

        if not Path(self.vcf).exists():
            raise ValueError(f"VCF文件不存在|VCF file does not exist: {self.vcf}")

        if not self.output:
            raise ValueError("输出VG文件不能为空|Output VG file cannot be empty")


@dataclass
class VGIndexConfig(VGBaseConfig):
    """VG Index配置类|VG Index Configuration Class"""

    # 输入文件|Input files
    input_graph: str = ""  # 输入图文件|Input graph file
    output_prefix: str = ""  # 输出前缀|Output prefix

    # 索引类型|Index types
    xg: bool = True  # 创建XG索引|Create XG index
    gcsa: bool = False  # 创建GCSA索引|Create GCSA index
    gbwt: bool = False  # 创建GBWT索引|Create GBWT index
    giraffe: bool = False  # 创建GIRAFFE索引|Create GIRAFFE indexes

    # GCSA参数|GCSA parameters
    kmer_size: int = 16  # k-mer大小|k-mer size
    edge_max: int = 4  # 最大边缘数|Maximum edges

    # 性能参数|Performance parameters
    threads: int = 12  # 线程数|Number of threads
    progress: bool = True  # 显示进度|Show progress

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        super().__post_init__()

        # 展开路径|Expand paths
        if self.input_graph:
            self.input_graph = expand_path(self.input_graph)
        if self.output_prefix:
            self.output_prefix = expand_path(self.output_prefix)

        # 确保输出目录存在|Ensure output directory exists
        if self.output_prefix:
            Path(self.output_prefix).parent.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        if not self.input_graph:
            raise ValueError("输入图文件不能为空|Input graph file cannot be empty")

        if not Path(self.input_graph).exists():
            raise ValueError(f"输入图文件不存在|Input graph file does not exist: {self.input_graph}")

        if not self.output_prefix:
            raise ValueError("输出前缀不能为空|Output prefix cannot be empty")


@dataclass
class VGGiraffeConfig(VGBaseConfig):
    """VG Giraffe配置类|VG Giraffe Configuration Class"""

    # 输入文件|Input files
    graph: str = ""  # 图文件前缀（索引）|Graph file prefix (indexed)
    reads: str = ""  # 输入reads文件|Input reads file
    output: str = ""  # 输出GAM文件|Output GAM file

    # 比对参数|Alignment parameters
    threads: int = 12  # 线程数|Number of threads
    fragment_length: int = 0  # 片段长度（0=自动检测）|Fragment length (0=auto)
    fragment_std_dev: int = 0  # 片段长度标准差（0=自动检测）|Fragment length std dev (0=auto)
    min_identity: float = 0.0  # 最小相似度|Min identity

    # 双端测序|Paired-end reads
    reads2: Optional[str] = None  # 第二个reads文件|Second reads file

    # 输出格式|Output format
    output_format: str = "GAM"  # 输出格式: GAM, GAF|Output format

    # 日志|Logging
    progress: bool = False  # 显示进度|Show progress

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        super().__post_init__()

        # 展开路径|Expand paths
        if self.graph:
            self.graph = expand_path(self.graph)
        if self.reads:
            self.reads = expand_path(self.reads)
        if self.reads2:
            self.reads2 = expand_path(self.reads2)
        if self.output:
            self.output = expand_path(self.output)

        # 规范化输出格式|Normalize output format
        self.output_format = self.output_format.upper()

        # 确保输出目录存在|Ensure output directory exists
        if self.output:
            Path(self.output).parent.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        if not self.graph:
            raise ValueError("图文件不能为空|Graph file cannot be empty")

        if not Path(self.graph).exists():
            raise ValueError(f"图文件不存在|Graph file does not exist: {self.graph}")

        if not self.reads:
            raise ValueError("Reads文件不能为空|Reads file cannot be empty")

        if not Path(self.reads).exists():
            raise ValueError(f"Reads文件不存在|Reads file does not exist: {self.reads}")

        if not self.output:
            raise ValueError("输出文件不能为空|Output file cannot be empty")


@dataclass
class VGDeconstructConfig(VGBaseConfig):
    """VG Deconstruct配置类|VG Deconstruct Configuration Class"""

    # 输入文件|Input files
    input_graph: str = ""  # 输入图文件|Input graph file
    output: str = ""  # 输出VCF文件|Output VCF file

    # 参考路径|Reference path
    reference_path: str = ""  # 参考路径名称|Reference path name

    # 样本选择|Sample selection
    samples: Optional[List[str]] = None  # 样本列表|Sample list

    # 输出参数|Output parameters
    threads: int = 12  # 线程数|Number of threads

    # 日志|Logging
    progress: bool = True  # 显示进度|Show progress

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        super().__post_init__()

        # 展开路径|Expand paths
        if self.input_graph:
            self.input_graph = expand_path(self.input_graph)
        if self.output:
            self.output = expand_path(self.output)

        # 确保输出目录存在|Ensure output directory exists
        if self.output:
            Path(self.output).parent.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        if not self.input_graph:
            raise ValueError("输入图文件不能为空|Input graph file cannot be empty")

        if not Path(self.input_graph).exists():
            raise ValueError(f"输入图文件不存在|Input graph file does not exist: {self.input_graph}")

        if not self.output:
            raise ValueError("输出VCF文件不能为空|Output VCF file cannot be empty")

        if not self.reference_path:
            raise ValueError("参考路径名称不能为空|Reference path name cannot be empty")
