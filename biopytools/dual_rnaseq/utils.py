"""
互作转录组工具函数模块|Dual RNA-seq Utility Functions Module
"""

import os
import logging
import subprocess
import sys
import time
from pathlib import Path


class DualRNASeqLogger:
    """互作转录组日志管理器|Dual RNA-seq Logger Manager"""

    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 创建日志文件|Create log file
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = output_dir / f"dual_rnaseq_processing_{timestamp}.log"

        # 配置logger|Configure logger
        self.logger = logging.getLogger(f"dual_rnaseq_processing_{timestamp}")

        # 设置日志级别|Set log level
        self.logger.setLevel(logging.DEBUG)

        # 清除现有的处理器|Clear existing handlers
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件处理器 - 记录所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        # stdout handler - INFO 及以下|stdout handler - INFO and below
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING 及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

    def get_logger(self):
        """获取logger实例|Get logger instance"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger):
        self.logger = logger

    def run(self, cmd: str, description: str = "", timeout: int = None) -> bool:
        """执行命令|Execute command

        Args:
            cmd: 要执行的命令|Command to execute
            description: 命令描述|Command description
            timeout: 超时时间（秒），None表示无限制|Timeout in seconds, None means no limit

        Returns:
            bool: 执行成功返回True，失败返回False|True if successful, False otherwise
        """
        if description:
            self.logger.info(f"运行|Running: {description}")

        self.logger.info(f"命令|Command: {cmd}")

        if timeout:
            self.logger.info(f"超时设置|Timeout: {timeout}秒|seconds ({timeout/3600:.1f}小时|hours)")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                timeout=timeout
            )

            self.logger.info(f"{description} 完成|completed")
            return True

        except subprocess.TimeoutExpired:
            self.logger.error(f"{description} 超时|timed out after {timeout}秒|seconds")
            return False

        except subprocess.CalledProcessError as e:
            self.logger.error(f"{description} 失败|failed")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False


class FileValidator:
    """文件验证器|File Validator"""

    def __init__(self, logger):
        self.logger = logger

    def check_file_exists(self, file_path: str, description: str = "") -> bool:
        """检查文件是否存在|Check if file exists

        Args:
            file_path: 文件路径|File path
            description: 文件描述|File description

        Returns:
            bool: 文件存在返回True，否则返回False|True if exists, False otherwise
        """
        if os.path.exists(file_path):
            if description:
                self.logger.info(f"{description}已存在，跳过|already exists, skipping: {file_path}")
            return True
        return False
