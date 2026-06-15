"""
KMC配置管理模块|KMC Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union
from ..common.paths import expand_path


@dataclass
class KMCConfig:
    """KMC配置类|KMC Configuration Class"""

    # 输入文件|Input files（支持单文件或目录|Support single file or directory）
    input_files: Optional[List[str]] = None
    input_dir: Optional[str] = None
    sample_names: Optional[List[str]] = None

    # 核心参数|Core parameters
    kmer_size: int = 21
    min_count: int = 2
    max_count: Optional[int] = None

    # 路径配置|Path configuration
    kmc_path: str = '~/miniforge3/envs/kmc_v.3.2.4/bin'
    output_dir: str = './kmc_output'
    tmp_dir: str = './kmc_tmp'

    # 处理参数|Processing parameters
    threads: int = 12
    memory_limit: Optional[str] = None  # 例如: "12G", "16G"

    # 双末端测序参数|Paired-end sequencing parameters
    read1_suffix: str = '_1.clean.fq.gz'  # read1文件后缀|read1 file suffix
    read2_suffix: str = '_2.clean.fq.gz'  # read2文件后缀|read2 file suffix
    single_end: bool = False  # 是否为单末端测序|Is single-end sequencing

    # 矩阵参数|Matrix parameters
    matrix_format: str = 'hdf5'  # 'hdf5', 'tsv', 'sqlite'
    sparse_storage: bool = True
    keep_dump: bool = True  # 是否保留dump文件|Whether to keep dump files
    max_memory: int = 500  # 最大内存使用量(GB)|Maximum memory usage (GB)

    # k-mer索引参数|K-mer index parameters
    index_mode: str = 'auto'  # 'auto', 'memory', 'db'
    index_threshold_gb: float = 1.0  # 自动选择阈值（GB）|Auto selection threshold (GB)

    # 操作模式|Operation mode
    mode: str = 'count'  # 'count', 'matrix', 'query', 'add'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开工具路径|Expand tool paths
        self.kmc_path = expand_path(self.kmc_path)

        # 标准化路径|Normalize paths（将相对路径转换为绝对路径|Convert relative paths to absolute paths）
        self.output_path = Path(self.output_dir).resolve()
        self.output_path.mkdir(parents=True, exist_ok=True)

        self.kmc_db_path = self.output_path / 'kmc_databases'
        self.kmc_db_path.mkdir(parents=True, exist_ok=True)

        self.tmp_path = Path(self.tmp_dir).resolve()
        self.tmp_path.mkdir(parents=True, exist_ok=True)

        # 创建 dump 文件目录|Create dump files directory
        self.dump_path = self.output_path / 'dump_files'
        self.dump_path.mkdir(parents=True, exist_ok=True)

        # 设置输入路径|Set input path
        if self.input_dir:
            self.input_path = Path(self.input_dir).resolve()
            if not self.input_path.exists():
                raise ValueError(f"输入目录不存在|Input directory does not exist: {self.input_dir}")
        elif self.input_files:
            self.input_path = None
        else:
            # 对于 matrix 和 query 模式，可以没有输入
            if self.mode not in ['matrix', 'query']:
                raise ValueError("必须指定 input_files 或 input_dir|Must specify either input_files or input_dir")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查KMC路径|Check KMC path
        if not os.path.exists(self.kmc_path):
            errors.append(f"KMC路径不存在|KMC path does not exist: {self.kmc_path}")

        kmc_bin = os.path.join(self.kmc_path, 'kmc')
        if not os.path.exists(kmc_bin):
            errors.append(f"KMC可执行文件不存在|KMC executable does not exist: {kmc_bin}")

        # 检查输入文件|Check input files（仅在指定了input_files时）
        if self.input_files:
            for input_file in self.input_files:
                if not os.path.exists(input_file):
                    errors.append(f"输入文件不存在|Input file does not exist: {input_file}")

        # 验证k-mer大小|Validate k-mer size
        if self.kmer_size <= 0 or self.kmer_size > 256:
            errors.append(f"k-mer大小必须在1-256之间|k-mer size must be between 1-256: {self.kmer_size}")

        # 验证最小计数|Validate minimum count
        if self.min_count < 1:
            errors.append(f"最小计数必须>=1|Minimum count must be >= 1: {self.min_count}")

        # 验证最大计数|Validate maximum count
        if self.max_count is not None and self.max_count < self.min_count:
            errors.append(
                f"最大计数({self.max_count})必须大于最小计数({self.min_count})|"
                f"Max count ({self.max_count}) must be greater than min count ({self.min_count})"
            )

        # 验证线程数|Validate thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        # 验证矩阵格式|Validate matrix format
        valid_formats = ['hdf5', 'tsv', 'sqlite']
        if self.matrix_format not in valid_formats:
            errors.append(
                f"无效的矩阵格式|Invalid matrix format: {self.matrix_format} "
                f"(必须是|must be one of {valid_formats})"
            )

        # 验证操作模式|Validate operation mode
        valid_modes = ['count', 'matrix', 'query', 'add']
        if self.mode not in valid_modes:
            errors.append(
                f"无效的操作模式|Invalid operation mode: {self.mode} "
                f"(必须是|must be one of {valid_modes})"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_kmc_bin(self):
        """获取KMC可执行文件路径|Get KMC executable path"""
        return os.path.join(self.kmc_path, 'kmc')

    def get_kmc_tools_bin(self):
        """获取kmc_tools可执行文件路径|Get kmc_tools executable path"""
        return os.path.join(self.kmc_path, 'kmc_tools')

    def get_kmc_dump_bin(self):
        """获取kmc_dump可执行文件路径|Get kmc_dump executable path"""
        return os.path.join(self.kmc_path, 'kmc_dump')

    def get_sample_db_path(self, sample_name: str) -> str:
        """获取样本数据库路径|Get sample database path"""
        return str(self.kmc_db_path / sample_name)

    def get_global_db_path(self) -> str:
        """获取全局数据库路径|Get global database path"""
        return str(self.output_path / 'global_kmers')
