"""Hi-C热图分析工具函数|Hi-C heatmap analysis utility functions"""

import logging
import sys
import subprocess
import os
import re
import shutil
from pathlib import Path
from ..common.conda_runner import build_conda_command  # §13 同源conda绝对路径+run -p, 严禁裸调conda
from typing import Optional, List


class HiCLogger:
    """Hi-C日志管理器|Hi-C Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers,
            force=True  # 确保重新配置|Ensure reconfiguration
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def check_dependencies(config, logger):
    """检查依赖工具|Check dependency tools

    Args:
        config: HiCConfig配置对象|HiCConfig object
        logger: 日志器|Logger

    Returns:
        bool: 是否所有工具都可用|Whether all tools are available
    """
    logger.info("检查依赖工具|Checking dependency tools...")

    tools = {
        'Juicer': config.juicer_sh,
        'Java': config.java_path,
        'BWA': config.bwa_path,
        'samtools': config.samtools_path,
        'plothic': config.plothic_path
    }

    missing = []

    for tool_name, tool_path in tools.items():
        try:
            # 针对不同工具使用不同的版本检查方法|Use different version check methods for different tools
            if tool_name == 'BWA':
                # BWA不支持--version，输出在stderr|BWA doesn't support --version, output in stderr
                result = subprocess.run(
                    [tool_path],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                # 检查stderr中是否包含版本信息|Check if stderr contains version info
                is_available = 'Version:' in result.stderr or 'Program: bwa' in result.stderr
            elif tool_name == 'Java':
                # Java使用-version参数|Java uses -version flag
                result = subprocess.run(
                    [tool_path, '-version'],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                is_available = 'version' in result.stderr.lower() or 'version' in result.stdout.lower()
            elif tool_name == 'Juicer':
                # Juicer脚本检查|Juicer script check
                is_available = os.path.exists(tool_path) and os.access(tool_path, os.X_OK)
            elif tool_name == 'plothic':
                # plothic命令检查|plothic command check
                result = subprocess.run(
                    [tool_path, '--version'],
                    capture_output=True,
                    text=True,
                    timeout=30  # 增加超时时间到30秒|Increase timeout to 30s
                )
                is_available = result.returncode == 0
            else:
                # 其他工具使用--version|Other tools use --version
                result = subprocess.run(
                    [tool_path, '--version'],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                is_available = result.returncode == 0

            if is_available:
                logger.debug(f"工具可用|Tool available: {tool_name}")
            else:
                missing.append(tool_name)
                logger.warning(f"工具未找到|Tool not found: {tool_name} at {tool_path}")

        except FileNotFoundError:
            missing.append(tool_name)
            logger.warning(f"工具未找到|Tool not found: {tool_name} at {tool_path}")
        except subprocess.TimeoutExpired:
            logger.warning(f"工具响应超时|Tool timeout: {tool_name}")
            missing.append(tool_name)
        except Exception as e:
            logger.warning(f"工具检查失败|Tool check failed: {tool_name} - {e}")
            missing.append(tool_name)

    if missing:
        logger.error(f"缺少依赖工具|Missing required tools: {', '.join(missing)}")
        return False

    logger.info("所有依赖工具可用|All dependency tools available")
    return True


def run_command(cmd, logger, description=""):
    """运行命令并记录日志|Run command and log

    Args:
        cmd: 命令列表|Command list
        logger: 日志器|Logger
        description: 命令描述|Command description

    Returns:
        bool: 是否成功|Whether successful
    """
    try:
        logger.debug(f"运行|Running: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )

        if description:
            logger.info(f"{description}完成|{description} completed")

        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败|Command execution failed: {e}")
        if e.stderr:
            logger.error(f"错误信息|Error message: {e.stderr}")
        return False

    except Exception as e:
        logger.error(f"命令执行异常|Command execution error: {e}")
        return False


def format_number(num: int) -> str:
    """格式化数字|Format number

    Args:
        num: 数字|Number

    Returns:
        str: 格式化后的字符串|Formatted string
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def get_file_size(file_path: str) -> str:
    """获取文件大小|Get file size

    Args:
        file_path: 文件路径|File path

    Returns:
        str: 格式化的文件大小|Formatted file size
    """
    try:
        size = os.path.getsize(file_path)
        return format_number(size)
    except Exception:
        return "Unknown"
