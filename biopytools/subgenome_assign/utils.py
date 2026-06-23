"""
亚基因组归属工具函数模块|Subgenome Assignment Utility Functions Module
"""

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple


class SubgenomeLogger:
    """亚基因组归属日志管理器|Subgenome Assignment Logger Manager"""

    def __init__(self, log_file: Path, log_name: str = "subgenome_assign.log"):
        self.log_dir = Path(log_file).parent
        self.log_file = self.log_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        formatter = logging.Formatter(log_format, datefmt=date_format)

        logger = logging.getLogger('biopytools.subgenome_assign')
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        # 文件handler|File handler (DEBUG+)
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        # stdout handler|Stdout handler (INFO+)
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler|Stderr handler (WARNING+)
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Optional[Path] = None):
        self.logger = logger
        self.working_dir = str(working_dir) if working_dir else None

    def run(self, cmd: List[str], description: str = "") -> bool:
        """执行单条命令（无管道）|Execute single command (no pipe)"""
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd, shell=False, capture_output=True, text=True,
                check=False, cwd=self.working_dir,
            )
            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command failed (rc={result.returncode}): {description}")
                if result.stderr:
                    self.logger.error(f"错误输出|Stderr: {result.stderr.strip()[:1000]}")
                return False
            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout.strip()[:500]}")
            self.logger.info(f"命令执行成功|Command succeeded: {description}")
            return True
        except Exception as e:
            self.logger.error(f"命令执行异常|Command exception: {description}")
            self.logger.error(f"异常信息|Exception: {e}")
            return False

    def run_pipeline(self, cmd_str: str, description: str = "") -> bool:
        """执行含重定向的命令|Execute command with redirection"""
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {cmd_str}")

        try:
            result = subprocess.run(
                cmd_str, shell=True, executable='/bin/bash',
                capture_output=True, text=True, check=False,
                cwd=self.working_dir,
            )
            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command failed (rc={result.returncode}): {description}")
                if result.stderr:
                    self.logger.error(f"错误输出|Stderr: {result.stderr.strip()[:1000]}")
                return False
            self.logger.info(f"命令执行成功|Command succeeded: {description}")
            return True
        except Exception as e:
            self.logger.error(f"命令执行异常|Command exception: {description}")
            self.logger.error(f"异常信息|Exception: {e}")
            return False


def check_dependencies(config, logger) -> bool:
    """检查依赖软件|Check dependencies"""
    logger.info("检查依赖软件|Checking dependencies")

    deps = [
        (config.minimap2_path, 'minimap2', ['--version']),
        (config.samtools_path, 'samtools', ['--version']),
    ]

    for path, name, vargs in deps:
        try:
            result = subprocess.run(
                [path] + vargs, capture_output=True, text=True, timeout=30
            )
            first_line = (result.stdout or result.stderr).strip().split('\n')[0]
            logger.info(f"{name} 可用|{name} available: {first_line}")
        except Exception as e:
            logger.error(f"{name} 不可用|{name} not available: {e}")
            return False
    return True


def format_number(n: int) -> str:
    """格式化大数字|Format large number"""
    if n >= 1_000_000:
        return f"{n / 1_000_000:.2f}M"
    if n >= 1_000:
        return f"{n / 1_000:.2f}K"
    return str(n)
