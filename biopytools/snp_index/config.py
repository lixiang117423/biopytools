"""
SNP Index配置管理模块 | SNP Index Configuration Management Module
"""

import os
import logging
from typing import Optional, Dict, Any


class SNPIndexConfig:
    """SNP Index分析配置类 | SNP Index Analysis Configuration Class"""

    def __init__(self,
                 input_vcf: Optional[str] = None,
                 output_file: Optional[str] = None,
                 output_dir: str = "./snp_index_output",
                 prefix: str = "snp_index",
                 min_depth: int = 10,
                 min_quality: int = 20,
                 min_mapping_quality: int = 20,
                 extreme_threshold: float = 0.8,
                 region_threshold: float = 0.5,
                 min_region_snps: int = 5,
                 max_region_gap: int = 10000,
                 sample_names: Optional[list] = None,
                 threads: int = 1,
                 keep_intermediate: bool = False,
                 force: bool = False,
                 quiet: bool = False,
                 verbose: int = 0,
                 log_file: Optional[str] = None,
                 # 滑动窗口参数 | Sliding window parameters
                 window_size: int = 1000000,
                 step_size: int = 100000,
                 min_window_snps: int = 5,
                 confidence_level: float = 0.95,
                 enable_sliding_window_plot: bool = True,
                 create_multi_chrom_plot: bool = False):
        """
        初始化配置 | Initialize configuration

        Args:
            input_vcf: 输入VCF文件路径 | Input VCF file path
            output_file: 输出文件路径 | Output file path
            output_dir: 输出目录 | Output directory
            prefix: 输出文件前缀 | Output file prefix
            min_depth: 最小测序深度 | Minimum sequencing depth
            min_quality: 最小质量值 | Minimum quality value
            min_mapping_quality: 最小mapping质量 | Minimum mapping quality
            extreme_threshold: 极端ΔSNP index阈值 | Extreme ΔSNP index threshold
            region_threshold: 区域检测阈值 | Region detection threshold
            min_region_snps: 区域最少SNP数量 | Minimum SNPs for region
            max_region_gap: 区域最大gap | Maximum gap in region
            sample_names: 样本名称列表 | Sample names list
            threads: 线程数 | Number of threads
            keep_intermediate: 保留中间文件 | Keep intermediate files
            force: 强制覆盖 | Force overwrite
            quiet: 静默模式 | Quiet mode
            verbose: 详细输出级别 | Verbose level
            log_file: 日志文件 | Log file
            window_size: 滑动窗口大小(bp) | Sliding window size in bp (default: 1000000)
            step_size: 滑动步长(bp) | Sliding step size in bp (default: 100000)
            min_window_snps: 窗口最少SNP数 | Minimum SNPs per window (default: 5)
            confidence_level: 置信水平 | Confidence level (default: 0.95)
            enable_sliding_window_plot: 启用滑动窗口折线图 | Enable sliding window line plot (default: True)
            create_multi_chrom_plot: 创建多染色体分离图 | Create multi-chromosome separated plot (default: False)
        """
        self.input_vcf = input_vcf
        self.output_file = output_file
        self.output_dir = output_dir
        self.prefix = prefix
        self.min_depth = min_depth
        self.min_quality = min_quality
        self.min_mapping_quality = min_mapping_quality
        self.extreme_threshold = extreme_threshold
        self.region_threshold = region_threshold
        self.min_region_snps = min_region_snps
        self.max_region_gap = max_region_gap
        self.sample_names = sample_names or []
        self.threads = threads
        self.keep_intermediate = keep_intermediate
        self.force = force
        self.quiet = quiet
        self.verbose = verbose
        self.log_file = log_file

        # 滑动窗口参数 | Sliding window parameters
        self.window_size = window_size
        self.step_size = step_size
        self.min_window_snps = min_window_snps
        self.confidence_level = confidence_level
        self.enable_sliding_window_plot = enable_sliding_window_plot
        self.create_multi_chrom_plot = create_multi_chrom_plot

        # 内部属性 | Internal attributes
        self._logger = None
        self._sample_columns = {}

    def validate(self) -> bool:
        """
        验证配置参数 | Validate configuration parameters

        Returns:
            bool: 验证是否通过 | Whether validation passed
        """
        if self.input_vcf and not os.path.exists(self.input_vcf):
            raise FileNotFoundError(f"输入VCF文件不存在 | Input VCF file not found: {self.input_vcf}")

        if self.min_depth < 0:
            raise ValueError("最小深度不能小于0 | Minimum depth cannot be negative")

        if not 0 <= self.extreme_threshold <= 1:
            raise ValueError("极端阈值必须在0-1之间 | Extreme threshold must be between 0-1")

        if not 0 <= self.region_threshold <= 1:
            raise ValueError("区域阈值必须在0-1之间 | Region threshold must be between 0-1")

        if self.threads < 1:
            raise ValueError("线程数必须大于0 | Threads must be greater than 0")

        return True

    def get_output_path(self, filename: str) -> str:
        """
        获取输出文件完整路径 | Get full output file path

        Args:
            filename: 文件名 | Filename

        Returns:
            str: 完整路径 | Full path
        """
        return os.path.join(self.output_dir, filename)

    def setup_output_dir(self) -> None:
        """创建输出目录 | Create output directory"""
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir, exist_ok=True)

    def get_log_level(self) -> int:
        """
        获取日志级别 | Get log level

        Returns:
            int: 日志级别 | Log level
        """
        if self.quiet:
            return logging.ERROR
        elif self.verbose >= 2:
            return logging.DEBUG
        elif self.verbose == 1:
            return logging.INFO
        else:
            return logging.WARNING

    def to_dict(self) -> Dict[str, Any]:
        """
        转换为字典 | Convert to dictionary

        Returns:
            dict: 配置字典 | Configuration dictionary
        """
        return {
            'input_vcf': self.input_vcf,
            'output_file': self.output_file,
            'output_dir': self.output_dir,
            'prefix': self.prefix,
            'min_depth': self.min_depth,
            'min_quality': self.min_quality,
            'min_mapping_quality': self.min_mapping_quality,
            'extreme_threshold': self.extreme_threshold,
            'region_threshold': self.region_threshold,
            'min_region_snps': self.min_region_snps,
            'max_region_gap': self.max_region_gap,
            'sample_names': self.sample_names,
            'threads': self.threads,
            'keep_intermediate': self.keep_intermediate,
            'force': self.force,
            'quiet': self.quiet,
            'verbose': self.verbose,
            'log_file': self.log_file
        }

    def __str__(self) -> str:
        """字符串表示 | String representation"""
        params = []
        for key, value in self.to_dict().items():
            params.append(f"{key}={value}")
        return f"SNPIndexConfig({', '.join(params)})"