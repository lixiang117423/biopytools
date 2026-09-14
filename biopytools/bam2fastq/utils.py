"""
BAM to FASTQ转换工具函数模块|BAM to FASTQ Conversion Utility Functions Module
"""

import logging
import sys
import subprocess
from pathlib import Path
from typing import Optional, List

# conda包装统一走公共层(§13): 同源conda绝对路径 + run -p <环境前缀>, 严禁裸调conda
from ..common.conda_runner import build_conda_command


class BAM2FASTQLogger:
    """BAM to FASTQ转换日志管理器|BAM to FASTQ Conversion Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "bam2fastq_conversion.log"):
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


class Bam2FQChecker:
    """bam2fastq工具检查器|bam2fastq Tool Checker"""

    def __init__(self, logger, bam2fastq_path: str = 'bam2fastq'):
        self.logger = logger
        self.bam2fastq_path = bam2fastq_path

    def check_bam2fastq(self):
        """检查bam2fastq是否已安装|Check if bam2fastq is installed"""
        try:
            cmd = build_conda_command(self.bam2fastq_path, ['--version'])
            self.logger.info(f"命令|Command: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            version = result.stdout.strip() if result.stdout.strip() else "unknown"
            self.logger.info(f"找到bam2fastq|Found bam2fastq: {version}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error(f"未找到bam2fastq|bam2fastq not found: {e}")
            self.logger.error("请先安装bam2fastq|Please install bam2fastq first")
            self.logger.error("安装方法|Installation: conda install -c bioconda pbmm2 pbbam")
            return False


def format_number(num: int) -> str:
    """格式化数字|Format number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)
