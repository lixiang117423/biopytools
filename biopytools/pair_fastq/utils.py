"""
FASTQ配对修复工具函数模块|FASTQ Pair Fixing Utility Functions Module
"""

import os
import sys
import logging
import subprocess
from pathlib import Path
from typing import Dict, Tuple, List
from datetime import datetime


class FastqPairLogger:
    """FASTQ配对修复日志管理器|FASTQ Pair Fixing Logger Manager"""

    def __init__(self, output_dir: Path, log_file: str = None, log_level: str = "INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_dir: 输出目录|Output directory
            log_file: 日志文件路径（可选）|Log file path (optional)
            log_level: 日志级别|Log level
        """
        self.output_dir = output_dir
        self.log_file = log_file
        self.log_level = getattr(logging, log_level.upper(), logging.INFO)

        # 创建日志目录|Create log directory
        self.log_dir = output_dir / "logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)

        # 设置默认日志文件|Set default log file
        if not self.log_file:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            self.log_file = self.log_dir / f"pair_fastq_{timestamp}.log"

        # 初始化logger|Initialize logger
        self.logger = self._setup_logger()

    def _setup_logger(self):
        """设置日志记录器|Setup logger"""
        logger = logging.getLogger("FastqPairFixer")
        logger.setLevel(self.log_level)

        # 清除已有的handlers|Clear existing handlers
        logger.handlers.clear()

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 控制台handler|Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(self.log_level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

        # 文件handler|File handler
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(self.log_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        return logger

    def get_logger(self):
        """获取logger对象|Get logger object"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志记录器|Logger object
        """
        self.logger = logger

    def run(self, command, description: str = None) -> bool:
        """
        执行命令|Execute command

        Args:
            command: 要执行的命令（字符串或列表）|Command to execute (string or list)
            description: 命令描述|Command description

        Returns:
            是否成功|Success status
        """
        if description:
            self.logger.info(f"执行|Executing: {description}")

        try:
            # 根据命令类型决定执行方式|Decide execution method based on command type
            if isinstance(command, list):
                # 列表形式：更安全，避免shell注入|List format: safer, no shell injection
                self.logger.debug(f"命令|Command: {' '.join(command)}")
                result = subprocess.run(
                    command,
                    check=False,  # 列表形式不使用shell，check=False以便捕获错误|Don't use shell for list
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
            else:
                # 字符串形式：使用shell执行|String format: execute with shell
                self.logger.debug(f"命令|Command: {command}")
                result = subprocess.run(
                    command,
                    shell=True,
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )

            # 检查返回码|Check return code
            if result.returncode != 0:
                self.logger.error(f"命令执行失败|Command execution failed: {description}")
                self.logger.error(f"返回码|Return code: {result.returncode}")
                if result.stderr:
                    self.logger.error(f"错误信息|Error: {result.stderr}")
                return False

            if result.stdout:
                self.logger.debug(f"输出|Output:\n{result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {description}")
            self.logger.error(f"返回码|Return code: {e.returncode}")
            if e.stderr:
                self.logger.error(f"错误信息|Error: {e.stderr}")
            return False

        except Exception as e:
            self.logger.error(f"执行命令时发生异常|Exception occurred while executing command: {str(e)}")
            return False


class FileManager:
    """文件管理器|File Manager"""

    @staticmethod
    def ensure_directory(dir_path: str):
        """
        确保目录存在|Ensure directory exists

        Args:
            dir_path: 目录路径|Directory path
        """
        Path(dir_path).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def find_files(directory: str, pattern: str) -> List[str]:
        """
        查找文件|Find files

        Args:
            directory: 搜索目录|Search directory
            pattern: 文件模式|File pattern

        Returns:
            文件路径列表|List of file paths
        """
        return list(Path(directory).glob(pattern))

    @staticmethod
    def get_file_size(file_path: str) -> str:
        """
        获取文件大小|Get file size

        Args:
            file_path: 文件路径|File path

        Returns:
            格式化的文件大小|Formatted file size
        """
        size_bytes = os.path.getsize(file_path)

        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if size_bytes < 1024.0:
                return f"{size_bytes:.2f} {unit}"
            size_bytes /= 1024.0

        return f"{size_bytes:.2f} PB"

    @staticmethod
    def count_files(directory: str, pattern: str = "*") -> int:
        """
        统计文件数量|Count files

        Args:
            directory: 目录路径|Directory path
            pattern: 文件模式|File pattern

        Returns:
            文件数量|Number of files
        """
        return len(list(Path(directory).glob(pattern)))


class PairFinder:
    """配对文件查找器|Paired Files Finder"""

    def __init__(self, input_dir: str, suffix1: str, suffix2: str, logger: logging.Logger):
        """
        初始化配对查找器|Initialize pair finder

        Args:
            input_dir: 输入目录|Input directory
            suffix1: R1文件后缀|R1 file suffix
            suffix2: R2文件后缀|R2 file suffix
            logger: 日志记录器|Logger
        """
        self.input_dir = input_dir
        self.suffix1 = suffix1
        self.suffix2 = suffix2
        self.logger = logger

    def find_pairs(self) -> Dict[str, Tuple[str, str]]:
        """
        查找配对文件|Find paired files

        Returns:
            字典: {sample_name: (r1_path, r2_path)}
            Dict: {sample_name: (r1_path, r2_path)}
        """
        self.logger.info(f"扫描输入目录: {self.input_dir}|Scanning input directory")

        r1_files = {}
        r2_files = {}

        # 查找R1文件|Find R1 files
        for file in os.listdir(self.input_dir):
            if file.endswith(self.suffix1):
                sample_name = file.replace(self.suffix1, '')
                r1_path = os.path.join(self.input_dir, file)
                r1_files[sample_name] = r1_path

        # 查找R2文件|Find R2 files
        for file in os.listdir(self.input_dir):
            if file.endswith(self.suffix2):
                sample_name = file.replace(self.suffix2, '')
                r2_path = os.path.join(self.input_dir, file)
                r2_files[sample_name] = r2_path

        # 找出完整配对|Find complete pairs
        pairs = {}
        r1_only = set(r1_files.keys()) - set(r2_files.keys())
        r2_only = set(r2_files.keys()) - set(r1_files.keys())

        for sample_name in r1_files.keys():
            if sample_name in r2_files:
                pairs[sample_name] = (r1_files[sample_name], r2_files[sample_name])

        # 报告统计信息|Report statistics
        self.logger.info(f"找到 {len(pairs)} 对完整文件|Found {len(pairs)} complete pairs")

        if r1_only:
            self.logger.warning(f"发现 {len(r1_only)} 个只有R1的样本|Found {len(r1_only)} R1-only samples")
            if len(r1_only) <= 5:
                self.logger.warning(f"  样本列表|Sample list: {list(r1_only)}")
            else:
                self.logger.warning(f"  部分样本|Partial samples: {list(r1_only)[:5]}...")

        if r2_only:
            self.logger.warning(f"发现 {len(r2_only)} 个只有R2的样本|Found {len(r2_only)} R2-only samples")
            if len(r2_only) <= 5:
                self.logger.warning(f"  样本列表|Sample list: {list(r2_only)}")
            else:
                self.logger.warning(f"  部分样本|Partial samples: {list(r2_only)[:5]}...")

        return pairs
