"""
JanusX工具模块|JanusX Utility Module

包含日志管理器和命令执行工具|Contains logger manager and command execution utilities
"""

import logging
import sys
import subprocess
from pathlib import Path
from typing import Optional, List


class JanusXLogger:
    """JanusX日志管理器|JanusX Logger Manager"""

    def __init__(self, output_path: Path, log_file: Optional[str] = None, log_level: str = "INFO"):
        """初始化日志管理器|Initialize logger manager

        Args:
            output_path: 输出目录路径|Output directory path
            log_file: 日志文件名|Log file name
            log_level: 日志级别|Log level
        """
        self.output_path = output_path
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level: str):
        """设置日志|Setup logging

        使用标准格式: YYYY-MM-DD HH:MM:SS.mmm - LEVEL - 中文|英文
        Use standard format: YYYY-MM-DD HH:MM:SS.mmm - LEVEL - Chinese|English
        """
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]

        # 添加文件处理器|Add file handler
        if self.log_file:
            log_path = self.output_path / self.log_file
            handlers.append(logging.FileHandler(log_path))

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers,
            force=True
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger, output_path: Path):
        """初始化命令执行器|Initialize command runner

        Args:
            logger: 日志器|Logger
            output_path: 输出目录路径|Output directory path
        """
        self.logger = logger
        self.output_path = output_path

    def run_command(self, cmd: List[str], description: str = "执行命令|Running command") -> bool:
        """运行命令|Run command

        Args:
            cmd: 命令列表|Command list
            description: 命令描述|Command description

        Returns:
            是否成功|Whether successful
        """
        self.logger.info(f"{description}: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )

            # 输出标准输出|Output stdout
            if result.stdout:
                for line in result.stdout.split('\n'):
                    if line.strip():
                        self.logger.info(f"输出|Output: {line}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {e}")

            # 输出错误信息|Output error message
            if e.stderr:
                for line in e.stderr.split('\n'):
                    if line.strip():
                        self.logger.error(f"错误|Error: {line}")

            return False

        except FileNotFoundError:
            self.logger.error(f"命令未找到|Command not found: {cmd[0]}")
            return False

        except Exception as e:
            self.logger.error(f"未知错误|Unknown error: {e}")
            return False


def check_janusx_dependencies(janusx_path: str, logger: logging.Logger):
    """检查JanusX依赖|Check JanusX dependencies

    Args:
        janusx_path: JanusX可执行文件路径|JanusX executable path
        logger: 日志器|Logger

    Raises:
        RuntimeError: 如果依赖检查失败|If dependency check fails
    """
    logger.info("检查JanusX依赖|Checking JanusX dependencies")

    errors = []

    # 检查JanusX可执行文件|Check JanusX executable
    if not janusx_path or not Path(janusx_path).exists():
        errors.append(f"JanusX可执行文件不存在|JanusX executable not found: {janusx_path}")
    else:
        logger.info(f"JanusX路径|JanusX path: {janusx_path}")

    # 检查JanusX是否可执行|Check if JanusX is executable
    if janusx_path and Path(janusx_path).exists():
        try:
            result = subprocess.run(
                [janusx_path, "-h"],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                logger.info("JanusX可执行文件检查通过|JanusX executable check passed")
            else:
                errors.append(f"JanusX执行失败|JanusX execution failed with code {result.returncode}")
        except subprocess.TimeoutExpired:
            errors.append("JanusX执行超时|JanusX execution timeout")
        except Exception as e:
            errors.append(f"JanusX执行错误|JanusX execution error: {e}")

    if errors:
        raise RuntimeError("\n".join(errors))
