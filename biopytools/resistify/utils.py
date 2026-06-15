"""
Resistify工具函数模块|Resistify Utility Functions Module
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name"""
    if os.path.isabs(command):
        match = re.search(r'/envs/([^/]+)', command)
        if match:
            return match.group(1)

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
    """构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment"""
    conda_env = get_conda_env(command)
    if conda_env:
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        full_cmd = [command] + args
    return full_cmd


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: str = None):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令|Execute command"""
        if description:
            self.logger.info(f"执行|Executing: {description}")

        self.logger.info(f"命令|Command: {cmd}")

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
            return False


class ResistifyLogger:
    """Resistify日志管理器|Resistify Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG/INFO/WARNING/ERROR)
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        level = getattr(logging, log_level.upper(), logging.INFO)

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler|File handler
        if self.log_file:
            file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger
