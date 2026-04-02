"""
Genome Analysis Configuration
基因组分析配置类
"""

import os


class GenomeAnalysisConfig:
    """基因组分析配置类|Genome Analysis Configuration Class"""

    def __init__(
        self,
        input_dir: str,
        output_dir: str,
        read_length: int = 150,
        kmer_size: int = 21,
        threads: int = 64,
        hash_size: str = "10G",
        max_kmer_cov: int = 1000,
        read1_suffix: str = "*_1.clean.fq.gz",
        skip_smudgeplot: bool = False,
        ploidy: int = 2,
        genomescope_env: str = "genomescope_v.2.0.1"
    ):
        """
        初始化配置|Initialize configuration

        Args:
            input_dir: 输入FASTQ文件目录|Input FASTQ directory
            output_dir: 输出目录|Output directory
            read_length: 测序读长 (默认: 150)|Read length (default: 150)
            kmer_size: K-mer大小 (默认: 21)|K-mer size (default: 21)
            threads: 线程数 (默认: 64)|Number of threads (default: 64)
            hash_size: Jellyfish哈希表大小 (默认: 10G)|Jellyfish hash size (default: 10G)
            max_kmer_cov: GenomeScope最大覆盖度 (默认: 1000)|Max k-mer coverage (default: 1000)
            read1_suffix: Read1文件后缀模式 (默认: *_1.clean.fq.gz)|Read1 file suffix pattern (default: *_1.clean.fq.gz)
            skip_smudgeplot: 跳过Smudgeplot倍性分析 (默认: False)|Skip Smudgeplot ploidy analysis (default: False)
            ploidy: 基因组倍性 1-6 (默认: 2)|Genome ploidy level 1-6 (default: 2)
            genomescope_env: GenomeScope conda环境名称 (默认: genomescope_v.2.0.1)|GenomeScope conda env name (default: genomescope_v.2.0.1)
        """
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.read_length = read_length
        self.kmer_size = kmer_size
        self.threads = threads
        self.hash_size = hash_size
        self.max_kmer_cov = max_kmer_cov
        self.read1_suffix = read1_suffix
        self.skip_smudgeplot = skip_smudgeplot
        self.ploidy = ploidy  # 用户指定的倍性（1-6），如果为None则由Smudgeplot推断
        self.genomescope_env = genomescope_env  # GenomeScope conda环境名称

        # 运行时变量|Runtime variables
        self.kcov = None  # GenomeScope计算得到的k-mer coverage
        self.inferred_ploidy = None  # Smudgeplot推断的倍性（如果未指定）

    def validate(self):
        """
        验证配置参数|Validate configuration parameters

        Raises:
            ValueError: 当配置参数无效时|When configuration parameters are invalid
        """
        # 检查必需目录|Check required directories
        if not self.input_dir or not os.path.exists(self.input_dir):
            raise ValueError(f"输入目录不存在|Input directory not found: {self.input_dir}")

        # 检查参数范围|Check parameter ranges
        if self.read_length <= 0:
            raise ValueError(f"读长必须大于0|Read length must be > 0: {self.read_length}")

        if self.kmer_size <= 0:
            raise ValueError(f"K-mer大小必须大于0|K-mer size must be > 0: {self.kmer_size}")

        if self.threads <= 0:
            raise ValueError(f"线程数必须大于0|Threads must be > 0: {self.threads}")

        if self.max_kmer_cov <= 0:
            raise ValueError(f"最大覆盖度必须大于0|Max coverage must be > 0: {self.max_kmer_cov}")

        if self.ploidy is not None and (self.ploidy < 1 or self.ploidy > 6):
            raise ValueError(f"倍性必须在1-6之间|Ploidy must be between 1-6: {self.ploidy}")

        return True

    def __repr__(self):
        """配置的字符串表示|String representation of configuration"""
        return (
            f"GenomeAnalysisConfig(\n"
            f"  input_dir={self.input_dir!r},\n"
            f"  output_dir={self.output_dir!r},\n"
            f"  read_length={self.read_length!r},\n"
            f"  kmer_size={self.kmer_size!r},\n"
            f"  threads={self.threads!r},\n"
            f"  hash_size={self.hash_size!r},\n"
            f"  max_kmer_cov={self.max_kmer_cov!r},\n"
            f"  read1_suffix={self.read1_suffix!r},\n"
            f"  skip_smudgeplot={self.skip_smudgeplot!r},\n"
            f"  ploidy={self.ploidy!r},\n"
            f"  genomescope_env={self.genomescope_env!r}\n"
            f")"
        )
