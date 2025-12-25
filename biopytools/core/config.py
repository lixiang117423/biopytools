"""
基础配置类 | Base Configuration Class
标准化配置基类，包含通用的参数验证和处理
"""

import os
import logging
from typing import Optional
from pathlib import Path


class BaseConfig:
    """基础配置类 | Base Configuration Class"""

    def __init__(self):
        """初始化基础配置"""
        self.version = "1.0.0"
        self.log_level = "INFO"
        self.verbose = 0
        self.quiet = False
        self.force = False
        self.dry_run = False
        self.keep_intermediate = False

    def validate_input_file(self, file_path: str, description: str = "输入文件") -> str:
        """
        验证输入文件存在性

        Args:
            file_path: 文件路径
            description: 文件描述

        Returns:
            验证后的文件路径

        Raises:
            FileNotFoundError: 文件不存在
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"{description}不存在 | {description} not found: {file_path}")
        return os.path.abspath(file_path)

    def validate_output_dir(self, output_dir: str, description: str = "输出目录") -> str:
        """
        验证并创建输出目录

        Args:
            output_dir: 输出目录路径
            description: 目录描述

        Returns:
            验证后的目录路径

        Raises:
            PermissionError: 无法创建目录
        """
        output_path = Path(output_dir)
        if not output_path.exists():
            output_path.mkdir(parents=True, exist_ok=True)
        elif not self.force and list(output_path.glob('*')):
            if not self.quiet:
                raise PermissionError(f"{description}非空，使用--force强制覆盖 | {description} not empty, use --force to overwrite")
        return str(output_path.absolute())

    def validate_threads(self, threads: int) -> int:
        """
        验证线程数

        Args:
            threads: 线程数

        Returns:
            验证后的线程数

        Raises:
            ValueError: 线程数无效
        """
        if threads < 1:
            raise ValueError(f"线程数必须大于0 | Thread count must be greater than 0: {threads}")
        return threads

    def validate_quality(self, quality: float, min_val: float = 0.0, max_val: float = 1.0) -> float:
        """
        验证质量参数

        Args:
            quality: 质量值
            min_val: 最小值
            max_val: 最大值

        Returns:
            验证后的质量值

        Raises:
            ValueError: 质量值超出范围
        """
        if not min_val <= quality <= max_val:
            raise ValueError(f"质量值必须在[{min_val}, {max_val}]范围内 | Quality value must be in range [{min_val}, {max_val}]: {quality}")
        return quality

    def get_log_level(self) -> int:
        """
        根据verbose和quiet设置获取日志级别

        Returns:
            logging level
        """
        if self.quiet:
            return logging.ERROR
        elif self.verbose >= 2:
            return logging.DEBUG
        elif self.verbose == 1:
            return logging.INFO
        else:
            return logging.WARNING

    def validate(self):
        """
        验证配置参数 - 子类应该重写此方法

        Raises:
            ValueError: 配置验证失败
        """
        pass

    def __str__(self) -> str:
        """配置的字符串表示 | String representation of configuration"""
        return f"""{self.__class__.__name__}:
    版本 | Version: {self.version}
    日志级别 | Log Level: {self.log_level}
    详细模式 | Verbose: {self.verbose}
    静默模式 | Quiet: {self.quiet}
    强制覆盖 | Force: {self.force}
    模拟运行 | Dry Run: {self.dry_run}
    保留中间文件 | Keep Intermediate: {self.keep_intermediate}
"""