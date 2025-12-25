"""
GTX Joint Calling工具类 | GTX Joint Calling Utility Classes
包含日志和辅助函数 | Contains logging and helper functions
"""

import logging
import sys
import os
import shutil
from pathlib import Path
from typing import List


class GTXJointLogger:
    """GTX Joint Calling日志管理器 | GTX Joint Calling Logger Manager"""

    def __init__(self, log_file: str = None):
        """
        初始化日志管理器 | Initialize logger manager

        Args:
            log_file: 日志文件路径 | Log file path
        """
        self.logger = logging.getLogger("GTXJoint")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()

        # 日志格式 | Log format
        formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO级别
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上级别
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler(如果指定) | File handler (if specified)
        if log_file:
            file_handler = logging.FileHandler(log_file, encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取logger对象 | Get logger object"""
        return self.logger


class GTXJointScanner:
    """GTX Joint Calling文件扫描器 | GTX Joint Calling File Scanner"""

    def __init__(self, logger):
        """
        初始化扫描器 | Initialize scanner

        Args:
            logger: 日志对象 | Logger object
        """
        self.logger = logger

    def scan_gvcf_files(self, gvcf_dir: str) -> List[str]:
        """
        扫描GVCF文件 | Scan GVCF files

        Args:
            gvcf_dir: GVCF目录 | GVCF directory

        Returns:
            GVCF文件列表 | List of GVCF files
        """
        self.logger.info(f"📂 正在扫描 GVCF 文件... | Scanning GVCF files...")
        gvcf_files = sorted(Path(gvcf_dir).glob("*.g.vcf.gz"))

        if not gvcf_files:
            self.logger.error(f"❌ 未找到任何 *.g.vcf.gz 文件 | No *.g.vcf.gz files found")
            return []

        # 检查索引文件
        missing_index = sum(1 for f in gvcf_files if not Path(f"{f}.tbi").exists())
        if missing_index > 0:
            self.logger.warning(
                f"⚠️ 有 {missing_index} 个文件缺少 .tbi 索引 | "
                f"{missing_index} files missing .tbi index"
            )

        self.logger.info(f"✅ 找到 {len(gvcf_files)} 个 GVCF 文件 | Found {len(gvcf_files)} GVCF files")
        return [str(f) for f in gvcf_files]

    def read_chromosomes(self, reference_fai: str) -> List[str]:
        """
        读取染色体列表 | Read chromosome list

        Args:
            reference_fai: 参考基因组索引文件 | Reference genome index file

        Returns:
            染色体列表 | List of chromosomes
        """
        self.logger.info(f"📖 正在读取染色体信息... | Reading chromosome information...")

        try:
            with open(reference_fai, 'r') as f:
                chromosomes = [line.split('\t')[0] for line in f]

            self.logger.info(f"✅ 参考基因组共有 {len(chromosomes)} 条染色体/scaffold | "
                           f"Reference has {len(chromosomes)} chromosomes/scaffolds")
            return chromosomes
        except Exception as e:
            self.logger.error(f"❌ 读取染色体信息失败 | Failed to read chromosome information: {e}")
            return []

    def check_faketime(self) -> bool:
        """
        检查faketime是否可用 | Check if faketime is available

        Returns:
            faketime是否可用 | Whether faketime is available
        """
        return shutil.which('faketime') is not None


class GTXJointWriter:
    """GTX Joint Calling命令写入器 | GTX Joint Calling Command Writer"""

    def __init__(self, logger):
        """
        初始化写入器 | Initialize writer

        Args:
            logger: 日志对象 | Logger object
        """
        self.logger = logger

    def backup_old_script(self, script_path: str) -> bool:
        """
        备份旧脚本 | Backup old script

        Args:
            script_path: 脚本路径 | Script path

        Returns:
            是否成功备份 | Whether backup was successful
        """
        if os.path.exists(script_path):
            backup = f"{script_path}.bak.{self._get_timestamp()}"
            try:
                shutil.move(script_path, backup)
                self.logger.warning(f"⚠️ 已备份旧文件为 | Old script backed up as: {backup}")
                return True
            except Exception as e:
                self.logger.error(f"❌ 备份失败 | Backup failed: {e}")
                return False
        return True

    def write_command(self, script_path: str, command: str) -> bool:
        """
        写入命令到脚本 | Write command to script

        Args:
            script_path: 脚本路径 | Script path
            command: 命令字符串 | Command string

        Returns:
            是否成功写入 | Whether write was successful
        """
        try:
            with open(script_path, 'a') as f:
                f.write(command + '\n')
            return True
        except Exception as e:
            self.logger.error(f"❌ 写入命令失败 | Failed to write command: {e}")
            return False

    def make_executable(self, script_path: str) -> bool:
        """
        设置脚本为可执行 | Make script executable

        Args:
            script_path: 脚本路径 | Script path

        Returns:
            是否成功设置 | Whether setting was successful
        """
        try:
            os.chmod(script_path, 0o755)
            return True
        except Exception as e:
            self.logger.error(f"❌ 设置可执行权限失败 | Failed to set executable permission: {e}")
            return False

    def _get_timestamp(self) -> str:
        """获取时间戳 | Get timestamp"""
        from datetime import datetime
        return datetime.now().strftime("%Y%m%d_%H%M%S")
