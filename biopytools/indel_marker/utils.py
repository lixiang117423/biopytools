"""INDEL分子标记工具函数模块|INDEL Marker Utility Functions Module"""

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional, List


class IndelMarkerLogger:
    """INDEL分子标记日志管理器|INDEL Marker Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "indel_marker.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        formatter = logging.Formatter(log_format, datefmt=date_format)

        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

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

    def _log_cmd(self, cmd: List[str], description: str):
        """记录命令到INFO级别|Log command at INFO level"""
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

    def run(self, cmd: List[str], description: str = "") -> bool:
        """
        执行命令（不捕获stdout，适合输出到文件的命令）

        Run command without returning stdout (suitable for file-output commands)

        Args:
            cmd: 命令列表（由build_conda_command()构建）|Command list (built by build_conda_command())
            description: 步骤描述|Step description

        Returns:
            bool: 执行是否成功|Whether execution succeeded
        """
        self._log_cmd(cmd, description)
        try:
            result = subprocess.run(cmd, shell=False, capture_output=True,
                                    text=True, check=False, cwd=self.working_dir)
            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command failed: {description}")
                if result.stderr:
                    self.logger.error(f"错误信息|Error: {result.stderr}")
                return False
            self.logger.info(f"命令执行成功|Command succeeded: {description}")
            return True
        except Exception as e:
            self.logger.error(f"命令执行异常|Command exception: {description}: {e}")
            return False

    def run_capture(self, cmd: List[str], description: str = "") -> Optional[str]:
        """
        执行命令并返回stdout|Run and return captured stdout

        Args:
            cmd: 命令列表（由build_conda_command()构建）|Command list (built by build_conda_command())
            description: 步骤描述|Step description

        Returns:
            stdout字符串，失败返回None|stdout string, or None on failure
        """
        self._log_cmd(cmd, description)
        try:
            result = subprocess.run(cmd, shell=False, capture_output=True,
                                    text=True, check=False, cwd=self.working_dir)
            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command failed: {description}")
                if result.stderr:
                    self.logger.error(f"错误信息|Error: {result.stderr}")
                return None
            return result.stdout
        except Exception as e:
            self.logger.error(f"命令执行异常|Command exception: {description}: {e}")
            return None


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令的conda环境|Detect conda env of a command

    策略|Strategy:
    1. 从which命令路径检测|Detect from which-path
    2. 兜底：搜索所有conda环境|Fallback: search all conda envs

    Args:
        command: 命令名称或完整路径|Command name or full path

    Returns:
        conda环境名称或None|Conda env name or None
    """
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    conda_exe = os.environ.get('CONDA_EXE')
    if conda_exe:
        conda_base_dir = os.path.dirname(os.path.dirname(conda_exe))
        envs_dir = os.path.join(conda_base_dir, 'envs')
        if os.path.exists(envs_dir):
            command_name = os.path.basename(command)
            for env_name in os.listdir(envs_dir):
                if os.path.exists(os.path.join(envs_dir, env_name, 'bin', command_name)):
                    return env_name
    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令|Build conda run command

    必须使用--no-capture-output避免conda缓冲输出导致内存问题
    Must use --no-capture-output to avoid conda buffering output causing memory issues

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        # 传完整路径(规范推荐,与phobius/braker统一)|Pass full path (recommended)
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


def check_dependencies(config, logger) -> bool:
    """
    检查依赖软件|Check dependencies（bcftools / samtools / primer3-py）

    Args:
        config: 配置对象，含bcftools_path/samtools_path|Config with bcftools_path/samtools_path
        logger: 日志器|Logger

    Returns:
        全部可用返回True，否则False|True if all available, False otherwise
    """
    logger.info("检查依赖软件|Checking dependencies")
    ok = True

    # primer3-py（Python包）|primer3-py Python package
    try:
        import primer3  # noqa: F401
        logger.info("primer3-py 可用|primer3-py available")
    except ImportError:
        logger.error("primer3-py 未安装，请 conda/pip install primer3-py|primer3-py not installed")
        ok = False

    # bcftools / samtools|external tools
    runner = CommandRunner(logger, Path.cwd())
    for tool_attr, label in [('bcftools_path', 'bcftools'), ('samtools_path', 'samtools')]:
        path = getattr(config, tool_attr)
        cmd = build_conda_command(path, ['--version'])
        out = runner.run_capture(cmd, description=f"{label} --version")
        if out is None:
            logger.error(f"{label} 不可用|{label} not available: {path}")
            ok = False
        else:
            logger.info(f"{label} 可用|{label} available: {out.strip().splitlines()[0] if out else ''}")
    return ok
