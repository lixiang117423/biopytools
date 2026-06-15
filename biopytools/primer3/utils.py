"""
Primer3引物设计工具函数模块|Primer3 Primer Design Utility Functions Module
"""

import logging
import sys
import subprocess
import re
import os
from pathlib import Path
from typing import Optional


class Primer3Logger:
    """Primer3引物设计日志管理器|Primer3 Primer Design Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "primer3_design.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
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


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    策略|Strategy:
    1. 首先尝试从which命令路径检测（优先级高）|First try detecting from which command path (high priority)
    2. 如果未找到，搜索所有conda环境（兜底方案）|If not found, search all conda environments (fallback)

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|conda environment name or None
    """
    import shutil

    # 方法1: 从命令路径检测|Method 1: Detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains envs
        # 例如: /miniforge3/envs/primer3_v.2.6.1/bin/primer3_core
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 方法2: 搜索所有conda环境|Method 2: Search all conda environments
    conda_base = os.environ.get('CONDA_EXE')
    if conda_base:
        # CONDA_EXE通常是/path/to/miniforge3/bin/conda|CONDA_EXE is usually /path/to/miniforge3/bin/conda
        conda_base_dir = os.path.dirname(os.path.dirname(conda_base))
        envs_dir = os.path.join(conda_base_dir, 'envs')

        if os.path.exists(envs_dir):
            # 搜索所有环境中的命令|Search command in all environments
            for env_name in os.listdir(envs_dir):
                env_bin = os.path.join(envs_dir, env_name, 'bin', command)
                if os.path.exists(env_bin):
                    return env_name

    return None


def build_conda_command(command: str, args: list) -> list:
    """
    构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Full command list
    """
    conda_env = get_conda_env(command)

    if conda_env:
        # 使用conda run调用|Use conda run to invoke
        full_cmd = ['conda', 'run', '-n', conda_env, command] + args
    else:
        # 非conda环境，直接调用|Non-conda environment, call directly
        full_cmd = [command] + args

    return full_cmd


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: list, description: str = "") -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command()构建）|Command list (built by build_conda_command())
            description: 步骤描述|Step description
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
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
            return False


def format_sequence_id(seq_id: str) -> str:
    """
    格式化序列ID，确保可以作为Primer3输入|Format sequence ID for Primer3 input

    Primer3要求SEQUENCE_ID不能包含特殊字符和空格|Primer3 requires SEQUENCE_ID without special chars and spaces

    Args:
        seq_id: 原始序列ID|Original sequence ID

    Returns:
        格式化后的序列ID|Formatted sequence ID
    """
    # 替换空格和特殊字符为下划线|Replace spaces and special chars with underscore
    formatted = re.sub(r'[^\w\-.]', '_', seq_id)
    return formatted


def parse_product_size_range(range_str: str) -> tuple:
    """
    解析产物大小范围字符串|Parse product size range string

    Args:
        range_str: 范围字符串，如"100-300"|Range string, e.g., "100-300"

    Returns:
        (min_size, max_size) 元组|(min_size, max_size) tuple
    """
    match = re.match(r'(\d+)-(\d+)', range_str)
    if match:
        return (int(match.group(1)), int(match.group(2)))
    return (100, 300)  # 默认值|Default
