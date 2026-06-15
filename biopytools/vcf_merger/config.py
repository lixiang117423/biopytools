"""
VCF合并工具配置模块|VCF Merger Configuration Module
"""

import os
from pathlib import Path
from typing import Optional


class VCFMergerConfig:
    """VCF合并工具配置类|VCF Merger Configuration Class"""

    def __init__(
        self,
        input_dir: str,
        output_dir: str,
        pattern: str = "*.joint.vcf.gz",
        threads: int = 4,
        create_index: bool = True,
        log_file: Optional[str] = None,
        log_level: str = "INFO",
        quiet: bool = False,
        verbose: int = 0
    ):
        """
        初始化配置|Initialize configuration

        Args:
            input_dir: 输入VCF文件目录|Input directory containing VCF files
            output_dir: 输出目录|Output directory
            pattern: VCF文件名模式|VCF file name pattern (default: *.joint.vcf.gz)
            threads: 线程数|Number of threads (default: 4)
            create_index: 是否创建索引文件|Whether to create index files (default: True)
            log_file: 日志文件路径|Log file path (optional)
            log_level: 日志级别|Log level (default: INFO)
            quiet: 静默模式|Quiet mode - only ERROR (default: False)
            verbose: 详细输出级别|Verbose level (0=WARNING, 1=INFO, 2=DEBUG)
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.pattern = pattern
        self.threads = threads
        self.create_index = create_index
        self.log_file = log_file
        self.log_level = log_level
        self.quiet = quiet
        self.verbose = verbose

    def validate(self):
        """
        验证配置参数|Validate configuration parameters

        Raises:
            ValueError: 配置参数无效时|When configuration parameters are invalid
            FileNotFoundError: 输入目录不存在时|When input directory does not exist
        """
        # 验证输入目录|Validate input directory
        if not self.input_dir.exists():
            raise FileNotFoundError(
                f"输入目录不存在|Input directory does not exist: {self.input_dir}"
            )

        if not self.input_dir.is_dir():
            raise ValueError(
                f"输入路径不是目录|Input path is not a directory: {self.input_dir}"
            )

        # 验证线程数|Validate thread count
        if self.threads < 1:
            raise ValueError(
                f"线程数必须大于0|Thread count must be greater than 0: {self.threads}"
            )

        # 验证日志级别|Validate log level
        valid_log_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        if self.log_level.upper() not in valid_log_levels:
            raise ValueError(
                f"无效的日志级别|Invalid log level: {self.log_level}. "
                f"必须是|Must be one of {valid_log_levels}"
            )

    def get_effective_log_level(self):
        """
        获取有效的日志级别|Get effective log level

        Returns:
            str: 有效的日志级别字符串|Effective log level string
        """
        if self.quiet:
            return "ERROR"
        elif self.verbose >= 2:
            return "DEBUG"
        elif self.verbose == 1:
            return "INFO"
        else:
            return self.log_level

    def __repr__(self):
        """配置的字符串表示|String representation of configuration"""
        return (
            f"VCFMergerConfig(\n"
            f"  input_dir={self.input_dir},\n"
            f"  output_dir={self.output_dir},\n"
            f"  pattern={self.pattern},\n"
            f"  threads={self.threads},\n"
            f"  create_index={self.create_index},\n"
            f"  log_level={self.get_effective_log_level()}\n"
            f")"
        )
