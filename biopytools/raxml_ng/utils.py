"""
RAxML-NG 分析工具函数模块|RAxML-NG Analysis Utility Functions Module
"""

import logging
import re
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


def get_conda_env(command: str) -> Optional[str]:
    """检测命令是否在conda环境中,返回环境名称|Detect if command is in a conda environment, return env name

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
    # 方法1: 从命令解析路径检测(绝对路径或 which 解析)|Method 1: detect from resolved path
    cmd_path = shutil.which(command)
    if cmd_path:
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    # 方法2: 仅当 command 是裸命令名(不含路径分隔符)时,才搜索所有 conda 环境|
    # Method 2: search all conda envs ONLY for bare command names (no path separator).
    # 对绝对路径(如静态二进制 /share/.../raxml-ng)必须跳过,否则 os.path.join 会塌缩到
    # 该绝对路径,误判为存在于某个 env 而返回错误的环境名。
    # For absolute paths (e.g. static binary) MUST skip, else os.path.join collapses
    # to that path and falsely matches the first env iterated.
    if os.path.sep not in command:
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

    静态二进制(不在 /envs/)直接调用;conda 环境内软件用 conda run --no-capture-output 包装
    Static binary (not under /envs/) called directly; conda env software wrapped with conda run --no-capture-output

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表(配合 subprocess.run(shell=False))|Complete command list (for subprocess.run(shell=False))
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        return [command] + args


class RaxmlNgLogger:
    """RAxML-NG 分析日志管理器|RAxML-NG Analysis Logger Manager"""

    def __init__(self, logs_dir: Path, log_name: str = "raxml_ng_analysis.log"):
        self.log_file = Path(logs_dir) / log_name
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self.setup_logging()

    def setup_logging(self):
        """设置日志(stdout INFO / stderr WARNING+ / 文件全级别)|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(log_format, datefmt=date_format)

        # stdout handler - INFO 级别 → 超算 .out|stdout handler - INFO → scheduler .out
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING 及以上 → 超算 .err|stderr handler - WARNING+ → scheduler .err
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件 handler - 所有级别 → 99_logs|file handler - all levels → 99_logs
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
        self.working_dir = Path(working_dir).resolve()

    def run(self, cmd: list, description: str = "") -> bool:
        """执行命令(列表,shell=False)|Execute command (list, shell=False)

        Args:
            cmd: 命令列表(由 build_conda_command 构建)|Command list (built by build_conda_command)
            description: 步骤描述|Step description

        Returns:
            成功 True,失败 False|True on success, False on failure
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        # 关键:记录完整命令到 INFO(便于论文 Methods 与可重复)|CRITICAL: log full command to INFO
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        self.logger.info(f"工作目录|Working directory: {self.working_dir}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir,
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


def check_dependencies(config, logger) -> Optional[str]:
    """检查依赖软件,返回版本字符串(失败返回None)|Check dependencies, return version string (None on failure)

    Args:
        config: RaxmlNgConfig 实例|RaxmlNgConfig instance
        logger: 日志器|Logger

    Returns:
        版本字符串或 None|Version string or None
    """
    logger.info("检查依赖软件|Checking dependencies")

    try:
        cmd = build_conda_command(config.raxml_ng_path, ['--version'])
        logger.info(f"命令|Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

        if result.returncode == 0:
            version_str = next(
                (ln.strip() for ln in (result.stdout or '').splitlines() if ln.strip()),
                'unknown',
            )
            logger.info(f"RAxML-NG 可用|RAxML-NG available: {version_str}")
            return version_str
        else:
            logger.error(f"RAxML-NG 版本检查失败|RAxML-NG version check failed (rc={result.returncode})")
            return None

    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        logger.error(f"缺少依赖软件|Missing dependency: RAxML-NG: {e}")
        logger.error("请确保 RAxML-NG 已安装或设置 RAXML_NG_PATH|Ensure RAxML-NG installed or set RAXML_NG_PATH")
        return None
