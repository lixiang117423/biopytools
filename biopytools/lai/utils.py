"""
LAI模块工具函数|LAI Module Utility Functions
"""

import logging
import sys
import subprocess
import os
import shutil
import re
from pathlib import Path
from typing import List, Optional


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 首先检查传入的command是否本身就是完整路径
    # First check if command itself is a full path
    if os.path.isabs(command):
        # 直接从完整路径中提取环境名|Extract env name directly from full path
        # 例如: /miniforge3/envs/ltr_harvest_parallel_v.1.2/bin/gt
        match = re.search(r'/envs/([^/]+)', command)
        if match:
            return match.group(1)

    # 如果不是完整路径，尝试用which查找|If not full path, try which
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains 'envs'
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 如果未找到，尝试搜索conda环境|If not found, try searching conda environments
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')

        if os.path.exists(envs_dir):
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment

    Args:
        command: 命令名称|Command name
        args: 命令参数|Command arguments

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, command] + args
    else:
        return [command] + args


class LAILogger:
    """LAI模块日志管理器|LAI Module Logger Manager"""

    def __init__(self, log_file: Path):
        self.log_file = log_file
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件handler|File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # stdout handler|Stdout handler
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[file_handler, stdout_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def run_command(cmd: list, logger: logging.Logger, check: bool = True,
                capture_output: bool = True) -> subprocess.CompletedProcess:
    """
    运行命令（自动检测conda环境）|Run command (auto-detect conda environment)

    Args:
        cmd: 命令列表|Command list
        logger: 日志器|Logger
        check: 是否检查返回码|Whether to check return code
        capture_output: 是否捕获输出|Whether to capture output

    Returns:
        subprocess.CompletedProcess: 命令执行结果|Command execution result
    """
    # 自动包装conda环境的命令|Auto-wrap conda environment commands
    if cmd:
        cmd_name = os.path.basename(cmd[0])
        wrapped_cmd = build_conda_command(cmd_name, cmd[1:])
    else:
        wrapped_cmd = cmd

    logger.debug(f"执行命令|Running command: {' '.join(wrapped_cmd)}")

    return subprocess.run(
        wrapped_cmd,
        capture_output=capture_output,
        text=True,
        check=check
    )


def check_dependencies(config, logger) -> bool:
    """
    检查依赖软件|Check dependencies

    根据运行模式检查所需软件|Check required tools based on run mode

    Args:
        config: 配置对象|Configuration object
        logger: 日志器|Logger

    Returns:
        bool: 是否所有依赖都存在|Whether all dependencies exist
    """
    logger.info("检查依赖软件|Checking dependencies")

    # 根据模式确定需要的软件|Determine required tools based on mode
    mode_requirements = {
        'full': [
            ('GenomeTools (gt)', config.gt_path),
            ('LTR_FINDER_parallel', config.ltr_finder_path),
            ('LTR_retriever', config.ltr_retriever_path),
            ('LAI', config.lai_path),
        ],
        'harvest': [
            ('GenomeTools (gt)', config.gt_path),
            ('LTR_FINDER_parallel', config.ltr_finder_path),
        ],
        'retrieve': [
            ('LTR_retriever', config.ltr_retriever_path),
        ],
        'calculate': [
            ('LAI', config.lai_path),
        ],
    }

    # 获取当前模式需要的依赖|Get dependencies for current mode
    required_tools = mode_requirements.get(config.mode, [])

    all_exist = True
    for name, path in required_tools:
        if Path(path).exists():
            logger.info(f"找到|Found {name}: {path}")
        else:
            logger.error(f"未找到|Not found {name}: {path}")
            all_exist = False

    if not all_exist:
        logger.error(f"部分依赖软件未找到 (模式: {config.mode})|Some dependencies not found (mode: {config.mode})")
        return False

    logger.info("所有依赖软件检查通过|All dependencies check passed")
    return True


def format_number(num: float, decimal: int = 2) -> str:
    """
    格式化数字|Format number

    Args:
        num: 数字|Number
        decimal: 小数位数|Decimal places

    Returns:
        str: 格式化后的字符串|Formatted string
    """
    return f"{num:.{decimal}f}"


def get_file_size(file_path: Path) -> str:
    """
    获取文件大小|Get file size

    Args:
        file_path: 文件路径|File path

    Returns:
        str: 格式化的文件大小|Formatted file size
    """
    size = file_path.stat().st_size

    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return f"{size:.2f} {unit}"
        size /= 1024.0

    return f"{size:.2f} PB"
