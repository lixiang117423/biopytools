"""
Merqury QV计算配置管理模块|Merqury QV Calculation Configuration Management Module
"""

import os
from ..common.paths import expand_path
import glob
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


@dataclass
class MerquryQVConfig:
    """Merqury QV计算配置类|Merqury QV Calculation Configuration Class"""

    # 输入路径|Input paths
    fastq_dir: str
    genome_file: str

    # 输出配置|Output configuration
    output_dir: str = "./merqury_qv_output"

    # 分析参数|Analysis parameters
    kmer_size: Optional[int] = None  # None表示自动选择|None means auto-select
    threads: int = 24

    # 软件路径|Software paths
    conda_env: str = "~/miniforge3/envs/merqury_v.1.3/bin/"
    merqury_path: Optional[str] = None  # None表示使用conda环境|None means use conda env

    # 数据类型|Data type
    data_type: str = "auto"  # auto, illumina, hifi

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 规范化路径|Normalize paths
        self.conda_env = expand_path(self.conda_env)
        self.fastq_dir = os.path.normpath(os.path.abspath(self.fastq_dir))
        self.genome_file = os.path.normpath(os.path.abspath(self.genome_file))

        # 设置MERQURY路径|Set MERQURY path
        if self.merqury_path is None:
            # 如果未指定，尝试从conda环境查找|Try to find from conda env if not specified
            self.merqury_path = self._find_merqury()

        # 发现FASTQ文件|Discover FASTQ files
        self.fastq_files = self._discover_fastq_files()

    def _find_merqury(self) -> str:
        """查找Merqury安装路径|Find Merqury installation path"""
        # 常见的Merqury安装位置|Common Merqury installation locations
        possible_paths = [
            Path(self.conda_env) / "../share/merqury",
            Path("~/tmp/test_merquery/merqury-master"),
            Path(os.environ.get("MERQURY", "")),
        ]

        for path in possible_paths:
            if path.exists() and (path / "merqury.sh").exists():
                return str(path)

        raise ValueError(
            f"无法找到Merqury安装|Cannot find Merqury installation. "
            f"请设置MERQURY环境变量或在config中指定merqury_path|"
            f"Please set MERQURY environment variable or specify merqury_path in config"
        )

    def _discover_fastq_files(self) -> List[str]:
        """
        发现FASTQ文件|Discover FASTQ files

        Returns:
            list: FASTQ文件路径列表|List of FASTQ file paths
        """
        path = Path(self.fastq_dir)
        fastq_files = []

        # 检查是否为单个文件|Check if it's a single file
        if path.is_file():
            # 支持的扩展名|Supported extensions
            valid_extensions = ['.fastq.gz', '.fq.gz', '.fastq', '.fq',
                                '.FASTQ.GZ', '.FQ.GZ', '.FASTQ', '.FQ']
            if any(str(path).endswith(ext) for ext in valid_extensions):
                return [str(path)]
            else:
                raise ValueError(
                    f"文件格式不支持|Unsupported file format: {path}. "
                    f"支持|Support: .fastq, .fq, .fastq.gz, .fq.gz"
                )

        # 处理目录|Handle directory
        # 支持的扩展名|Supported extensions
        extensions = ['*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq']

        for ext in extensions:
            pattern = str(path / ext)
            fastq_files.extend(glob.glob(pattern))

        # 同时支持大写扩展名|Also support uppercase extensions
        for ext in ['*.FASTQ.GZ', '*.FQ.GZ', '*.FASTQ', '*.FQ']:
            pattern = str(path / ext)
            fastq_files.extend(glob.glob(pattern))

        return sorted(list(set(fastq_files)))  # 去重并排序|Remove duplicates and sort

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入路径|Check input paths
        if not os.path.exists(self.fastq_dir):
            errors.append(f"FASTQ目录不存在|FASTQ directory does not exist: {self.fastq_dir}")

        if not os.path.exists(self.genome_file):
            errors.append(f"基因组文件不存在|Genome file does not exist: {self.genome_file}")

        # 检查是否找到FASTQ文件|Check if FASTQ files were found
        if not self.fastq_files:
            errors.append(f"未找到FASTQ文件|No FASTQ files found in: {self.fastq_dir}")

        # 检查基因组格式|Check genome file format
        if not self.genome_file.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz')):
            errors.append(
                f"基因组文件格式不正确|Genome file format incorrect: {self.genome_file}. "
                f"支持.fa/.fasta格式|Support .fa/.fasta format"
            )

        # 检查参数范围|Check parameter ranges
        if self.kmer_size is not None and (self.kmer_size < 15 or self.kmer_size > 31):
            errors.append(f"Kmer大小必须在15-31之间|Kmer size must be between 15-31: {self.kmer_size}")

        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive: {self.threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
