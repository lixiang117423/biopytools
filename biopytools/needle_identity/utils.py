"""
NeedleIdentity工具函数模块|Needle Identity Utility Functions Module
"""

import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

# conda包装统一走公共层(§13): 同源conda绝对路径 + run -p <环境前缀>, 严禁裸调conda
from ..common.conda_runner import build_conda_command


class NeedleIdentityLogger:
    """NeedleIdentity日志管理器|Needle Identity Logger Manager"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level: str):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        self.logger = logging.getLogger('needle_identity')
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False

        formatter = logging.Formatter(log_format, datefmt=date_format)

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
            log_path = Path(self.log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def expand_path(path: str) -> str:
    """展开路径中的~和环境变量|Expand ~ and environment variables in path"""
    return os.path.expandvars(os.path.expanduser(path))


def get_tool_path(tool_name: str, default_path: str, env_var: Optional[str] = None) -> str:
    """
    按优先级获取工具路径|Get tool path by priority

    优先级(高到低)|Priority (high to low):
        1. 环境变量 (env_var)|Environment variable
        2. 用户配置文件 (~/.config/biopytools/config.yml)|User config file
        3. 代码默认值 (default_path)|Code default
    """
    # 1. 环境变量优先|Env var takes priority
    if env_var and os.environ.get(env_var):
        return expand_path(os.environ[env_var])

    # 2. 用户配置文件|User config file
    config_path = os.path.expanduser("~/.config/biopytools/config.yml")
    if os.path.exists(config_path):
        try:
            import yaml
            with open(config_path) as f:
                config = yaml.safe_load(f) or {}
            tool_path = config.get("tools", {}).get(tool_name)
            if tool_path:
                return expand_path(tool_path)
        except Exception:
            pass

    # 3. 代码默认值|Code default
    return expand_path(default_path)


def check_needle(config, logger: logging.Logger) -> bool:
    """检查EMBOSS needle可用性|Check EMBOSS needle availability"""
    logger.info("检查依赖软件|Checking dependencies: EMBOSS needle")
    embossversion_path = os.path.join(os.path.dirname(config.needle_path), 'embossversion')
    try:
        cmd = build_conda_command(embossversion_path, [])
        logger.info(f"命令|Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, shell=False, capture_output=True, text=True, timeout=60)
        if result.returncode == 0 and result.stdout.strip():
            version = result.stdout.strip().splitlines()[0]
            logger.info(f"EMBOSS版本|EMBOSS version: {version}")
            return True
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
        pass
    logger.error(f"EMBOSS不可用|EMBOSS not available: {config.needle_path}")
    return False
