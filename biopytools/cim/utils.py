"""
CIM分析工具函数模块|CIM Analysis Utility Functions Module
"""

import logging
import re
import shutil
import subprocess
import sys
import os
from typing import List, Optional, Tuple


class CIMLogger:
    """CIM分析日志管理器|CIM Analysis Logger Manager"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG/INFO/WARNING/ERROR)
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level: str):
        """设置日志|Setup logging"""
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
            force=True
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self) -> logging.Logger:
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志器|Logger instance
        """
        self.logger = logger

    def run(self, cmd: str, description: str = "",
            timeout: Optional[int] = None) -> Tuple[bool, str, str]:
        """
        执行命令|Execute command

        Args:
            cmd: 要执行的命令|Command to execute
            description: 命令描述|Command description
            timeout: 超时时间(秒)|Timeout in seconds

        Returns:
            tuple: (是否成功|success, 标准输出|stdout, 标准错误|stderr)
        """
        if description:
            self.logger.info(f"执行|Running: {description}")
        self.logger.info(f"命令|Command: {cmd}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=timeout
            )

            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command failed (exit {result.returncode})")
                if result.stderr:
                    # 只记录最后20行错误信息|Only log last 20 lines of stderr
                    err_lines = result.stderr.strip().split('\n')
                    for line in err_lines[-20:]:
                        self.logger.error(f"  {line}")
                return False, result.stdout, result.stderr

            return True, result.stdout, result.stderr

        except subprocess.TimeoutExpired:
            self.logger.error(f"命令执行超时|Command timed out ({timeout}s): {cmd}")
            return False, "", "Timeout"
        except Exception as e:
            self.logger.error(f"命令执行异常|Command error: {e}")
            return False, "", str(e)

    def run_conda(self, env_name: str, cmd: str, description: str = "",
                  timeout: Optional[int] = None) -> Tuple[bool, str, str]:
        """
        在conda环境中执行命令|Execute command in conda environment

        Args:
            env_name: conda环境名|Conda environment name
            cmd: 要执行的命令|Command to execute
            description: 命令描述|Command description
            timeout: 超时时间(秒)|Timeout in seconds

        Returns:
            tuple: (是否成功|success, 标准输出|stdout, 标准错误|stderr)
        """
        if '/' in env_name or '~' in env_name:
            full_cmd = f"conda run --no-capture-output -p {env_name} {cmd}"
        else:
            full_cmd = f"conda run --no-capture-output -n {env_name} {cmd}"
        return self.run(full_cmd, description, timeout)


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中,返回环境名称|Detect if command is in a conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

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
    """构建conda run命令来运行conda环境中的软件|Build conda run command for conda env software

    Args:
        command: 命令名称或完整路径(建议传完整路径以正确检测env)|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list

    Note:
        必须使用--no-capture-output避免conda缓冲输出导致大数据内存溢出
        |Must use --no-capture-output to avoid conda buffering causing OOM on large data
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args
