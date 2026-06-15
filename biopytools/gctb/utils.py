"""
GCTB工具函数模块|GCTB Utility Functions Module
"""

import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional, List


class GCTBLogger:
    """GCTB日志管理器|GCTB Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "gctb_pipeline.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 删除旧日志|Delete old log
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        # 文件handler|File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        # stdout handler|Stdout handler
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        # stderr handler|Stderr handler (WARNING+)
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        # 配置根日志|Configure root logger
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.DEBUG)
        root_logger.handlers.clear()
        root_logger.propagate = False

        root_logger.addHandler(file_handler)
        root_logger.addHandler(stdout_handler)
        root_logger.addHandler(stderr_handler)

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: List[str], description: str = "") -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command()构建）|Command list (built by build_conda_command())
            description: 步骤描述|Step description

        Returns:
            bool: 执行是否成功|Whether execution succeeded
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        # 记录完整命令到INFO级别|Log complete command at INFO level
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,  # 传入列表时必须使用shell=False|Must use shell=False with list
                capture_output=True,
                text=True,
                check=False,
                cwd=self.working_dir
            )

            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command execution failed: {description}")
                self.logger.error(f"错误代码|Error code: {result.returncode}")
                if result.stderr:
                    self.logger.error(f"错误信息|Error message: {result.stderr}")
                return False

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")
            return True

        except Exception as e:
            self.logger.error(f"命令执行异常|Command execution error: {description}")
            self.logger.error(f"异常信息|Exception: {str(e)}")
            return False


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
    import shutil
    import re

    # 方法1: 从命令路径检测|Method 1: Detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 方法2: 搜索所有conda环境|Method 2: Search all conda environments
    conda_exe = os.environ.get('CONDA_EXE')
    if conda_exe:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_exe))
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
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list

    ⚠️ 重要|IMPORTANT:
        必须使用--no-capture-output避免conda缓冲输出导致内存问题
        Must use --no-capture-output to avoid conda buffering output causing memory issues
    """
    import os

    # 从路径中提取命令名称|Extract command name from path
    command_name = os.path.basename(command)

    conda_env = get_conda_env(command)

    if conda_env:
        # 使用conda run调用，只使用命令名称|Use conda run with command name only
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command_name] + args
    else:
        # 非conda环境，使用完整路径或命令名称|Non-conda environment, use full path or command name
        full_cmd = [command] + args

    return full_cmd
