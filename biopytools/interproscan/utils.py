"""
InterProScan注释工具函数模块|InterProScan Annotation Utility Functions Module
"""

import logging
import os
import subprocess
import sys
from pathlib import Path


class InterProScanLogger:
    """InterProScan注释日志管理器|InterProScan Annotation Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "interproscan_annotation.log"):
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


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger, working_dir: Path = None, env_vars: dict = None):
        self.logger = logger
        self.working_dir = working_dir or Path.cwd()
        self.env_vars = env_vars or {}

    def run(self, cmd: str, description: str = "", env_vars: dict = None) -> bool:
        """执行命令|Execute command

        Args:
            cmd: 命令字符串|Command string
            description: 步骤描述|Step description
            env_vars: 环境变量字典|Environment variables dict (will be merged with init env_vars)
        """
        if description:
            self.logger.info(f"执行步骤|Executing step: {description}")

        self.logger.info(f"命令|Command: {cmd}")

        # 合并环境变量|Merge environment variables
        import os
        env = os.environ.copy()
        if self.env_vars:
            env.update(self.env_vars)
        if env_vars:
            env.update(env_vars)

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir,
                env=env
            )

            self.logger.info(f"命令执行成功|Command executed successfully: {description}")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False


def count_sequences(fasta_file: str) -> int:
    """统计FASTA文件序列数量|Count sequences in FASTA file"""
    count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def format_number(num: int) -> str:
    """格式化数字|Format number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)
