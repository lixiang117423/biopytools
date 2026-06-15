"""
K-mer计数配置管理模块|K-mer Count Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List


@dataclass
class KmerCountConfig:
    """K-mer计数配置类|K-mer Count Configuration Class"""

    # 输入文件|Input files
    input_dir: Path = field(default_factory=lambda: Path.cwd())
    pattern: str = "*_1.fq.gz"
    kmer_lib: Path = field(default_factory=lambda: Path.cwd() / "kmers.fasta")
    bed_file: Optional[Path] = None

    # 输出设置|Output settings
    output_dir: Path = field(default_factory=lambda: Path.cwd() / "kmer_count_output")

    # Jellyfish参数|Jellyfish parameters
    kmer_size: int = 51
    hash_size: str = '1000M'
    threads: int = 8
    canonical: bool = False

    # 滑动窗口参数|Sliding window parameters
    window_size: int = 500000
    step_size: Optional[int] = None

    # 其他参数|Other parameters
    keep_temp: bool = True
    keep_binary: bool = False
    verbose: bool = False

    # 工具路径|Tool paths
    jellyfish_path: str = 'jellyfish'

    # 内部变量|Internal variables
    temp_dir: Optional[Path] = field(default=None)
    samples: List[str] = field(default_factory=list)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 转换路径为Path对象|Convert paths to Path objects
        if isinstance(self.input_dir, str):
            self.input_dir = Path(self.input_dir)
        if isinstance(self.kmer_lib, str):
            self.kmer_lib = Path(self.kmer_lib)
        if isinstance(self.bed_file, str):
            self.bed_file = Path(self.bed_file)
        if isinstance(self.output_dir, str):
            self.output_dir = Path(self.output_dir)

        # 创建输出目录|Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 设置默认步长|Set default step size
        if self.step_size is None:
            self.step_size = self.window_size // 5

    def setup_temp_dir(self) -> Path:
        """设置临时目录|Setup temporary directory"""
        import time
        import random
        temp_name = f"kmer_count_{int(time.time())}_{random.randint(1000, 9999)}"
        self.temp_dir = Path.cwd() / "tmp" / temp_name
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        return self.temp_dir

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入目录|Check input directory
        if not self.input_dir.exists():
            errors.append(f"输入目录不存在|Input directory does not exist: {self.input_dir}")

        # 检查K-mer库文件|Check K-mer library file
        if not self.kmer_lib.exists():
            errors.append(f"K-mer库文件不存在|K-mer library file does not exist: {self.kmer_lib}")

        # 检查BED文件|Check BED file
        if self.bed_file and not self.bed_file.exists():
            errors.append(f"BED文件不存在|BED file does not exist: {self.bed_file}")

        # 检查参数范围|Check parameter ranges
        if self.kmer_size <= 0:
            errors.append(f"K-mer大小必须为正数|K-mer size must be positive: {self.kmer_size}")

        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if self.window_size <= 0:
            errors.append(f"窗口大小必须为正数|Window size must be positive: {self.window_size}")

        if self.step_size and self.step_size <= 0:
            errors.append(f"步长必须为正数|Step size must be positive: {self.step_size}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    @classmethod
    def from_args(cls, args) -> 'KmerCountConfig':
        """从命令行参数创建配置|Create configuration from command line arguments"""
        return cls(
            input_dir=args.input,
            pattern=args.pattern,
            kmer_lib=args.kmer_lib,
            bed_file=getattr(args, 'bed_file', None),
            output_dir=args.output,
            kmer_size=args.kmer_size,
            hash_size=args.hash_size,
            threads=args.threads,
            canonical=args.canonical,
            window_size=args.window_size,
            step_size=getattr(args, 'step_size', None),
            keep_temp=args.keep_temp,
            keep_binary=args.keep_binary,
            verbose=args.verbose,
            jellyfish_path=args.jellyfish_path
        )
