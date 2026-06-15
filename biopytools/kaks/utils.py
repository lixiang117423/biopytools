"""
Ka/Ks Calculator工具函数模块|Ka/Ks Calculator Utility Functions Module
"""

import logging
import re
import shutil
import sys
import os
from pathlib import Path
from typing import List, Optional


class KaKsLogger:
    """Ka/Ks分析日志管理器|Ka/Ks Analysis Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "kaks_analysis.log", verbose: bool = False):
        """
        初始化日志器|Initialize logger

        Args:
            output_dir: 输出目录(日志自动放到99_logs/子目录)|Output directory (log auto-placed in 99_logs/)
            log_name: 日志文件名|Log file name
            verbose: 详细模式|Verbose mode
        """
        self.output_dir = output_dir
        self.log_dir = output_dir / "99_logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.log_file = self.log_dir / log_name
        self._verbose = verbose
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        logger = logging.getLogger("KaKsAnalyzer")
        logger.setLevel(logging.DEBUG if self._verbose else logging.INFO)
        logger.handlers.clear()
        logger.propagate = False

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
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        self.logger = logger

    def info(self, message: str):
        """信息日志|Info logging"""
        self.logger.info(message)

    def success(self, message: str):
        """成功日志|Success logging"""
        self.logger.info(message)

    def warning(self, message: str):
        """警告日志|Warning logging"""
        self.logger.warning(message)

    def error(self, message: str):
        """错误日志|Error logging"""
        self.logger.error(message)

    def debug(self, message: str):
        """调试日志|Debug logging"""
        self.logger.debug(message)

    def progress(self, message: str, current: int, total: int):
        """进度日志|Progress logging"""
        percentage = (current / total) * 100 if total > 0 else 0
        self.logger.info(f"{message} [{current}/{total}] ({percentage:.1f}%)")

    def separator(self, title: str = ""):
        """分隔符日志|Separator logging"""
        if title:
            self.logger.info(f"{title} " + "="*50)
        else:
            self.logger.info("="*60)

    def get_logger(self):
        """获取日志器对象|Get logger object"""
        return self.logger


def get_conda_env(command: str) -> Optional[str]:
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda env, return env name

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
    构建conda run命令来运行conda环境中的软件|Build conda run command to invoke software in conda env

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Complete command list
    """
    conda_env = get_conda_env(command)
    if conda_env:
        return ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    return [command] + args


def format_number(num: int) -> str:
    """格式化大数字|Format large number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return f"{num:,}"
