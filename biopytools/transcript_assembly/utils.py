"""
转录本从头组装工具函数模块|Transcript De Novo Assembly Utility Functions Module
"""

import os
import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional, List

# conda包装统一走公共层(§13): 同源conda绝对路径 + run -p <环境前缀>, 严禁裸调conda
from ..common.conda_runner import build_conda_command


class TranscriptAssemblyLogger:
    """转录本组装日志管理器|Transcript Assembly Logger Manager"""

    def __init__(self, output_dir: Path, verbose: bool = False, quiet: bool = False,
                 log_name: str = "pipeline.log"):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 创建日志目录和文件|Create log directory and file
        log_dir = output_dir / "99_logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        self.log_file = log_dir / log_name

        # 配置logger|Configure logger
        self.logger = logging.getLogger("transcript_assembly")

        # 设置日志级别|Set log level
        if quiet:
            self.logger.setLevel(logging.ERROR)
        elif verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        # 清除现有的处理器|Clear existing handlers
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)

        self.logger.propagate = False

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取logger实例|Get logger instance"""
        return self.logger

    def step(self, message: str):
        """记录步骤分隔|Log step separator"""
        self.logger.info("=" * 60)
        self.logger.info(message)
        self.logger.info("=" * 60)


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger):
        self.logger = logger

    def run(self, cmd: str, description: str = "", timeout: int = None) -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 要执行的命令字符串|Command string to execute
            description: 命令描述|Command description
            timeout: 超时时间（秒）|Timeout in seconds
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        self.logger.info(f"命令|Command: {cmd}")

        if timeout:
            self.logger.info(f"超时设置|Timeout: {timeout}秒|seconds ({timeout / 3600:.1f}小时|hours)")

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
            self.logger.error(
                f"{description} 超时|timed out after {timeout}秒|seconds ({timeout / 3600:.1f}小时|hours)")
            self.logger.error(f"跳过该步骤继续处理|Skipping this step and continuing...")
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
        """检查文件是否存在|Check if file exists"""
        if os.path.exists(file_path):
            if description:
                self.logger.info(f"{description}已存在，跳过|already exists, skipping: {file_path}")
            return True
        return False
