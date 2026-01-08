"""
FASTP质控工具函数模块|FASTP Quality Control Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional


class FastpLogger:
    """FASTP质控日志管理器|FASTP Quality Control Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "fastp_processing.log",
                 log_level: str = "INFO", quiet: bool = False):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_dir: 输出目录路径|Output directory path
            log_name: 日志文件名|Log file name
            log_level: 日志级别 (DEBUG/INFO/WARNING/ERROR/CRITICAL)|Log level
            quiet: 静默模式（只输出ERROR）| Quiet mode (ERROR only)
        """
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.log_level = log_level
        self.quiet = quiet
        self.setup_logging()

    def setup_logging(self):
        """设置日志系统|Setup logging system"""
        # 确保输出目录存在|Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 创建logger|Create logger
        self.logger = logging.getLogger(__name__)
        self.logger.handlers.clear()
        self.logger.setLevel(logging.DEBUG)

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - 输出INFO及以下级别|stdout handler - INFO and below
        stdout_handler = logging.StreamHandler(sys.stdout)
        if self.quiet:
            # 静默模式：只输出ERROR到stderr，stdout不输出|Quiet mode: only ERROR to stderr
            stdout_handler.setLevel(logging.CRITICAL + 1)  # 禁用stdout
        else:
            stdout_handler.setLevel(logging.INFO)
            stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - 输出WARNING及以上级别|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler - 记录所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志对象|Logger object
        """
        self.logger = logger

    def run(self, cmd: list, description: str = "") -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表|Command list
            description: 步骤描述|Step description

        Returns:
            bool: 执行是否成功|Whether execution succeeded
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        # 将命令列表转换为字符串用于日志|Convert command list to string for logging
        cmd_str = " ".join(str(arg) for arg in cmd)
        self.logger.info(f"命令|Command: {cmd_str}")

        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            # 如果有输出，记录到debug级别|Log output to debug level if exists
            if result.stdout.strip():
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False

    def check_executable(self, executable_path: str, timeout: int = 10) -> bool:
        """
        检查可执行文件是否可用|Check if executable is available

        Args:
            executable_path: 可执行文件路径|Executable file path
            timeout: 超时时间（秒）| Timeout in seconds

        Returns:
            bool: 是否可用|Whether available
        """
        try:
            result = subprocess.run(
                [executable_path, "--version"],
                capture_output=True,
                text=True,
                timeout=timeout
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
            return False
