"""
Swave工具函数模块|Swave Utility Functions Module
"""

import logging
import sys
import subprocess
import os
import shutil
import re
from typing import Optional, List, Tuple


class SwaveLogger:
    """Swave日志管理器|Swave Logger Manager"""

    def __init__(self, output_dir: str, log_file: str = 'swave.log', log_level: str = "INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_dir: 输出目录|Output directory
            log_file: 日志文件名|Log file name
            log_level: 日志级别|Log level
        """
        self.log_file = log_file
        self.log_level = log_level
        self.output_dir = output_dir

        # 创建日志目录|Create log directory
        log_dir = os.path.join(output_dir, '99_logs')
        os.makedirs(log_dir, exist_ok=True)
        self.log_path = os.path.join(log_dir, log_file)

        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 清除现有的handlers|Clear existing handlers
        root_logger = logging.getLogger()
        root_logger.handlers.clear()

        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, self.log_level.upper(), logging.INFO)

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))

        # 配置root logger|Configure root logger
        root_logger.setLevel(level)
        root_logger.addHandler(stdout_handler)
        root_logger.addHandler(stderr_handler)
        root_logger.addHandler(file_handler)

        # 禁止传播|Disable propagation
        root_logger.propagate = False

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment, return env name

    策略|Strategy:
    1. 首先尝试从which命令路径检测（优先级高）|Try detecting from which command path first (high priority)
    2. 如果未找到，搜索所有conda环境（兜底方案）|If not found, search all conda environments (fallback)

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|conda environment name or None
    """
    # 方法1: 从命令路径检测|Method 1: Detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains envs
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


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件|Build conda run command to run software in conda environment

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list

    Examples:
        >>> build_conda_command('busco', ['--version'])
        ['conda', 'run', '-n', 'BUSCO_v.6.0.0', 'busco', '--version']

        >>> # 绝对路径且不在conda envs目录下时，直接调用
        >>> # Absolute path not under conda envs: called directly
        >>> build_conda_command('/usr/bin/tool', ['--help'])
        ['/usr/bin/tool', '--help']
    """
    conda_env = get_conda_env(command)

    if conda_env:
        # 使用conda run调用，--no-capture-output避免缓冲输出导致内存问题
        # Use conda run, --no-capture-output to avoid memory issues from output buffering
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        # 非conda环境，直接调用|Non-conda environment, call directly
        full_cmd = [command] + args

    return full_cmd


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger, output_dir: str):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志器|Logger
            output_dir: 输出目录|Output directory
        """
        self.logger = logger
        self.output_dir = output_dir

    def run(self, cmd: List[str], description: str = "",
            cwd: Optional[str] = None, env: Optional[dict] = None) -> Tuple[bool, str, str]:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表|Command list
            description: 步骤描述|Step description
            cwd: 工作目录|Working directory
            env: 环境变量|Environment variables

        Returns:
            (成功状态, 标准输出, 标准错误)|(Success status, stdout, stderr)
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        try:
            self.logger.debug(f"命令|Command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=False,
                cwd=cwd,
                env=env
            )

            if result.returncode != 0:
                self.logger.error(f"命令失败|Command failed: {description}")
                self.logger.error(f"错误输出|Error output: {result.stderr}")
                return False, result.stdout, result.stderr

            if description:
                self.logger.info(f"完成|Completed: {description}")

            return True, result.stdout, result.stderr

        except Exception as e:
            self.logger.error(f"执行异常|Execution error: {description}")
            self.logger.error(f"异常信息|Exception: {str(e)}")
            return False, "", str(e)


def check_dependencies(config, logger: logging.Logger) -> bool:
    """
    检查依赖软件|Check dependencies

    Args:
        config: Swave配置对象|Swave config object
        logger: 日志器|Logger

    Returns:
        是否所有依赖都可用|Whether all dependencies are available
    """
    logger.info("检查依赖软件|Checking dependencies")

    dependencies = []

    # 检查Python环境|Check Python environment
    swave_py = os.path.join(config.swave_path, 'Swave.py')
    if os.path.exists(swave_py):
        logger.info(f"Swave主程序|Swave main script: {swave_py}")
        dependencies.append(True)
    else:
        logger.error(f"Swave主程序不存在|Swave main script not found: {swave_py}")
        dependencies.append(False)

    # 检查minigraph（如果需要）|Check minigraph (if needed)
    if config.gfa_source == 'minigraph':
        minigraph_cmd = build_conda_command(config.minigraph_path, ['--version'])
        try:
            result = subprocess.run(minigraph_cmd, capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"minigraph可用|minigraph available: {config.minigraph_path}")
                dependencies.append(True)
            else:
                logger.error(f"minigraph不可用|minigraph not available: {config.minigraph_path}")
                dependencies.append(False)
        except Exception as e:
            logger.error(f"minigraph检测失败|minigraph check failed: {e}")
            dependencies.append(False)

    # 检查gfatools（可选）|Check gfatools (optional)
    try:
        gfatools_cmd = build_conda_command(config.gfatools_path, [])
        result = subprocess.run(gfatools_cmd, capture_output=True, text=True, timeout=10)
        # gfatools不支持--version，无参数时输出usage（含"gfatools"）并返回非0
        # gfatools doesn't support --version, outputs usage containing "gfatools" without args
        if result.returncode != 127 and 'gfatools' in result.stderr:
            logger.info(f"gfatools可用|gfatools available: {config.gfatools_path}")
        else:
            logger.warning(f"gfatools不可用（可选）|gfatools not available (optional)")
    except Exception:
        logger.warning(f"gfatools检测失败（可选）|gfatools check failed (optional)")

    all_ok = all(dependencies)
    if all_ok:
        logger.info("依赖检查完成|Dependency check completed")
    else:
        logger.error("部分依赖缺失|Some dependencies missing")

    return all_ok
