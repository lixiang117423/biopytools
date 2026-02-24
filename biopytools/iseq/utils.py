"""
iSeq下载工具函数模块|iSeq Download Utility Functions Module
"""

import logging
import os
import subprocess
import sys
from pathlib import Path


class ISeqLogger:
    """iSeq下载日志管理器|iSeq Download Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "iseq_download.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件handler|File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # stdout handler|Stdout handler
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[file_handler, stdout_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


class CondaCommandRunner:
    """支持conda环境的命令执行器|Command Runner with conda environment support"""

    def __init__(self, logger, conda_env: str, working_dir: Path):
        self.logger = logger
        self.conda_env = conda_env
        self.working_dir = working_dir

    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令（在conda环境中）|Execute command (in conda environment)"""
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {cmd}")

        # 使用conda run在指定环境中执行命令|Use conda run to execute in specified environment
        full_cmd = f"conda run -n {self.conda_env} --live-stream {cmd}"

        try:
            result = subprocess.run(
                full_cmd,
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
