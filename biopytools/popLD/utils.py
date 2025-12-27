"""
PopLDdecay工具函数模块 | PopLDdecay Utility Functions Module
"""

import logging
import os
import subprocess
import sys
from pathlib import Path
from datetime import datetime


class PopLDdecayLogger:
    """PopLDdecay日志管理器 | PopLDdecay Logger Manager"""

    def __init__(self, output_dir: str = ".", log_name: str = "poplddecay.log"):
        self.output_dir = Path(output_dir)
        self.log_file = self.output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志 | Setup logging"""
        # 创建日志目录 | Create log directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 配置日志
        log_format = '[%(asctime)s] %(levelname)s: %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        # 创建formatter
        formatter = logging.Formatter(log_format, datefmt=date_format)

        # 创建根logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()

        # stdout handler - INFO级别及以下
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上级别
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler - 所有级别
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器 | Command Runner"""

    def __init__(self, logger, working_dir: str = "."):
        self.logger = logger
        self.working_dir = working_dir

    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")

        self.logger.info(f"命令 | Command: {cmd}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功 | Command executed successfully: {description}")

            if result.stdout:
                # 只输出重要信息到DEBUG级别
                for line in result.stdout.strip().split('\n'):
                    if line.strip():
                        self.logger.debug(f"  {line}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            if e.stderr:
                for line in e.stderr.strip().split('\n'):
                    if line.strip():
                        self.logger.error(f"  {line}")
            return False

    def run_with_output(self, cmd: str, description: str = "") -> tuple:
        """执行命令并返回输出 | Execute command and return output"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")

        self.logger.info(f"命令 | Command: {cmd}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir
            )

            self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            return (True, result.stdout)

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            if e.stderr:
                for line in e.stderr.strip().split('\n'):
                    if line.strip():
                        self.logger.error(f"  {line}")
            return (False, e.stderr if e.stderr else str(e))
