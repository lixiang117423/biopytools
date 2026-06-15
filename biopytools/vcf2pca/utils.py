"""
VCF2PCA分析工具函数模块|VCF2PCA Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path


class VCF2PCALogger:
    """VCF2PCA日志管理器|VCF2PCA Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """
        设置日志|Setup logging

        Args:
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            # 确保日志目录存在|Ensure log directory exists
            log_path = Path(self.log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            handlers.append(logging.FileHandler(self.log_file))

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


# 保留旧类名以保持兼容性|Keep old class name for compatibility
class PCALogger(VCF2PCALogger):
    """PCA分析日志管理器（兼容性保留）|PCA Analysis Logger Manager (Legacy)"""

    def __init__(self, output_dir: Path, log_name: str = "pca_analysis.log"):
        """
        初始化（兼容旧接口）|Initialize (legacy interface)

        Args:
            output_dir: 输出目录|Output directory
            log_name: 日志文件名|Log file name
        """
        log_file = output_dir / log_name
        super().__init__(log_file=log_file, log_level="INFO")


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志器|Logger
            working_dir: 工作目录|Working directory
        """
        self.logger = logger
        self.working_dir = working_dir.resolve()  # 使用绝对路径|Use absolute path

    def run(self, cmd: str, description: str = "") -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 命令字符串|Command string
            description: 步骤描述|Step description

        Returns:
            bool: 是否成功|Success or not
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {cmd}")
        self.logger.info(f"工作目录|Working directory: {self.working_dir}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            self.logger.error(f"标准输出|Stdout: {e.stdout}")
            return False


def check_dependencies(config, logger):
    """
    检查依赖软件|Check dependencies

    Args:
        config: 配置对象|Configuration object
        logger: 日志器|Logger

    Returns:
        bool: 是否所有依赖都可用|All dependencies available

    Raises:
        RuntimeError: 缺少依赖时|Missing dependencies
    """
    logger.info("检查依赖软件|Checking dependencies")

    dependencies = [
        (config.plink_path, "PLINK"),
        (config.bcftools_path, "BCFtools")
    ]

    missing_deps = []

    for cmd, name in dependencies:
        try:
            result = subprocess.run(
                [cmd, "--version"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                logger.info(f"{name} 可用|{name} available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)

    if missing_deps:
        error_msg = f"缺少依赖软件|Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    return True
