"""
IQ-TREE分析工具函数模块|IQ-TREE Analysis Utility Functions Module
"""

import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中，返回环境名称|Detect if command is in a conda environment

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
    import re
    import os

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
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        return [command] + args


class TreeLogger:
    """系统发育树分析日志管理器|Phylogenetic Tree Analysis Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "iqtree_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(log_format, datefmt=date_format)

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()

    def run(self, cmd: list, description: str = "") -> bool:
        """执行命令|Execute command

        Args:
            cmd: 命令列表（由build_conda_command构建）|Command list (built by build_conda_command)
            description: 步骤描述|Step description
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        self.logger.info(f"工作目录|Working directory: {self.working_dir}")

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
            self.logger.error(f"标准输出|Stdout: {e.stdout}")
            return False


def check_dependencies(config, logger):
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    try:
        cmd = build_conda_command(config.iqtree_path, ['--version'])
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )

        if result.returncode == 0:
            version_line = result.stdout.split('\n')[0]
            logger.info(f"IQ-TREE 可用|IQ-TREE available: {version_line}")
            return True
        else:
            raise RuntimeError("IQ-TREE 未正确安装|IQ-TREE not properly installed")

    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        error_msg = f"缺少依赖软件|Missing dependency: IQ-TREE"
        logger.error(error_msg)
        logger.error(f"请确保IQ-TREE已安装并在PATH中|Please ensure IQ-TREE is installed and in PATH")
        raise RuntimeError(error_msg)
