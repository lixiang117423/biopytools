"""
ADMIXTURE分析工具函数模块|ADMIXTURE Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import shutil
from pathlib import Path
from typing import Optional, List


class AdmixtureLogger:
    """ADMIXTURE分析日志管理器|ADMIXTURE Analysis Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "admixture_analysis.log",
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

        # 清理旧的处理器，避免日志重复打印|Clear old handlers to avoid duplicate logging
        logger = logging.getLogger()
        if logger.hasHandlers():
            logger.handlers.clear()

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
        logger.addHandler(stdout_handler)

        # stderr handler - 输出WARNING及以上级别|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件handler - 记录所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        # 设置根日志级别|Set root logger level
        logger.setLevel(logging.DEBUG)
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志对象|Logger object
            working_dir: 工作目录|Working directory
        """
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: str, description: str = "") -> str:
        """
        执行命令|Execute command

        Args:
            cmd: 命令字符串|Command string
            description: 步骤描述|Step description

        Returns:
            命令输出|Command output
        """
        if description:
            self.logger.info(f"开始|Starting: {description}")

        cleaned_cmd = " ".join(cmd.strip().split())
        self.logger.info(f"执行命令|Executing command: {cleaned_cmd}")

        try:
            result = subprocess.run(
                cleaned_cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                cwd=self.working_dir
            )

            if result.stdout.strip():
                self.logger.info(f"命令输出|Command output:\n{result.stdout.strip()}")

            if result.stderr.strip():
                self.logger.warning(f"命令警告|Command warnings:\n{result.stderr.strip()}")

            self.logger.info(f"完成|Completed: {description}")
            return result.stdout

        except subprocess.CalledProcessError as e:
            self.logger.error(
                f"命令执行失败|Command execution failed: {cleaned_cmd}\n"
                f"   - 返回码|Return code: {e.returncode}\n"
                f"   - 标准输出|Stdout: {e.stdout.strip()}\n"
                f"   - 标准错误|Stderr: {e.stderr.strip()}"
            )
            raise


class SoftwareChecker:
    """软件环境检查器|Software Environment Checker"""

    def __init__(self, logger):
        """
        初始化软件检查器|Initialize software checker

        Args:
            logger: 日志对象|Logger object
        """
        self.logger = logger

    def check_software(self, software: str) -> bool:
        """
        检查软件是否安装|Check if software is installed

        Args:
            software: 软件名称|Software name

        Returns:
            是否已安装|Whether installed
        """
        if shutil.which(software):
            self.logger.info(f"{software} 已安装|{software} is installed")
            return True
        else:
            self.logger.warning(f"{software} 未安装|{software} is not installed")
            return False

    def check_dependencies(self) -> bool:
        """
        检查所有依赖软件|Check all dependencies

        Returns:
            是否所有必需软件都已安装|Whether all required software installed
        """
        required_software = ["plink", "admixture", "bcftools"]
        optional_software = ["Rscript"]

        self.logger.info("=== 软件环境检查|Software Environment Check ===")

        all_required = True
        for software in required_software:
            if not self.check_software(software):
                all_required = False

        for software in optional_software:
            self.check_software(software)

        if not all_required:
            self.logger.error("缺少必需软件，请安装后再运行|Missing required software, please install and try again")
            return False

        return True


def check_file_format(file_path: str, expected_format: str) -> bool:
    """
    检查文件格式|Check file format

    Args:
        file_path: 文件路径|File path
        expected_format: 期望格式 (vcf/bed)|Expected format

    Returns:
        格式是否匹配|Whether format matches
    """
    if expected_format == "vcf":
        return file_path.endswith(('.vcf', '.vcf.gz'))
    elif expected_format == "bed":
        return file_path.endswith('.bed')
    return False
