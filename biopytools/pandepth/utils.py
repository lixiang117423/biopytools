"""
PanDepth覆盖度计算工具函数模块|PanDepth Coverage Calculation Utility Functions Module
"""

import os
import logging
import subprocess
import sys
import time
from pathlib import Path


class PanDepthLogger:
    """PanDepth覆盖度计算日志管理器|PanDepth Coverage Calculation Logger Manager"""

    def __init__(self, output_dir: Path, verbose: bool = False, quiet: bool = False, log_name: str = "pandepth_processing.log"):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 创建日志文件|Create log file
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = output_dir / f"pandepth_processing_{timestamp}.log"

        # 配置logger|Configure logger
        self.logger = logging.getLogger(f"pandepth_processing_{timestamp}")

        # 设置日志级别|Set log level
        if quiet:
            self.logger.setLevel(logging.ERROR)
        elif verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        # 清除现有的处理器|Clear existing handlers
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件处理器 - 记录所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        # stdout handler - INFO 及以下|stdout handler - INFO and below
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING 及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

    def get_logger(self):
        """获取logger实例|Get logger instance"""
        return self.logger

    def step(self, message: str):
        """记录步骤|Log step"""
        self.logger.info("=" * 60)
        self.logger.info(message)
        self.logger.info("=" * 60)


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger):
        self.logger = logger

    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令|Execute command

        Args:
            cmd: 要执行的命令|Command to execute
            description: 命令描述|Command description
        """
        if description:
            self.logger.info(f"运行|Running: {description}")

        self.logger.info(f"命令|Command: {cmd}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )

            self.logger.info(f"{description} 完成|completed")

            if result.stdout:
                self.logger.debug(f"标准输出|Stdout: {result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"{description} 失败|failed")
            self.logger.error(f"错误代码|Error code: {e.returncode}")
            self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False


class BAMFileFinder:
    """BAM文件查找器|BAM File Finder"""

    def __init__(self, logger):
        self.logger = logger

    def find_bam_files(self, input_path: str) -> list:
        """查找BAM文件|Find BAM files

        Args:
            input_path: 输入路径(文件或目录)|Input path (file or directory)

        Returns:
            list: BAM文件路径列表|List of BAM file paths
        """
        bam_files = []

        if os.path.isdir(input_path):
            # 目录模式：查找所有BAM文件|Directory mode: find all BAM files
            self.logger.info(f"扫描目录查找BAM文件|Scanning directory for BAM files: {input_path}")

            for file in Path(input_path).glob('*.bam'):
                # 跳过BAM索引文件和其他扩展名|.bam.bai, .bam.csi等
                if not file.name.endswith('.bam'):
                    continue
                bam_files.append(str(file))

            if not bam_files:
                self.logger.warning(f"目录中未找到BAM文件|No BAM files found in directory: {input_path}")

            self.logger.info(f"找到 {len(bam_files)} 个BAM文件|Found {len(bam_files)} BAM files")

        elif os.path.isfile(input_path):
            # 单文件模式|Single file mode
            if input_path.endswith('.bam'):
                bam_files.append(input_path)
                self.logger.info(f"使用单个BAM文件|Using single BAM file: {input_path}")
            else:
                self.logger.error(f"输入文件不是BAM格式|Input file is not in BAM format: {input_path}")

        return bam_files
