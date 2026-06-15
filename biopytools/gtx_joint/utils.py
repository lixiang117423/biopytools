"""
GTX Joint Calling工具类|GTX Joint Calling Utility Classes
包含日志和辅助函数|Contains logging and helper functions
"""

import logging
import sys
import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional


class GTXJointLogger:
    """GTX Joint Calling日志管理器|GTX Joint Calling Logger Manager"""

    def __init__(self, log_file: str = None):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
        """
        self.logger = logging.getLogger("GTXJoint")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
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

        # 文件handler(如果指定)|File handler (if specified)
        if log_file:
            file_handler = logging.FileHandler(log_file, encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取logger对象|Get logger object"""
        return self.logger


class GTXJointScanner:
    """GTX Joint Calling文件扫描器|GTX Joint Calling File Scanner"""

    def __init__(self, logger):
        """
        初始化扫描器|Initialize scanner

        Args:
            logger: 日志对象|Logger object
        """
        self.logger = logger

    def scan_gvcf_files(self, gvcf_dir: str) -> List[str]:
        """
        扫描GVCF文件|Scan GVCF files

        Args:
            gvcf_dir: GVCF目录|GVCF directory

        Returns:
            GVCF文件列表|List of GVCF files
        """
        self.logger.info(f"正在扫描 GVCF 文件...|Scanning GVCF files...")
        gvcf_files = sorted(Path(gvcf_dir).glob("*.g.vcf.gz"))

        if not gvcf_files:
            self.logger.error(f"未找到任何 *.g.vcf.gz 文件|No *.g.vcf.gz files found")
            return []

        # 检查索引文件
        missing_index = sum(1 for f in gvcf_files if not Path(f"{f}.tbi").exists())
        if missing_index > 0:
            self.logger.warning(
                f"有 {missing_index} 个文件缺少 .tbi 索引|"
                f"{missing_index} files missing .tbi index"
            )

        self.logger.info(f"找到 {len(gvcf_files)} 个 GVCF 文件|Found {len(gvcf_files)} GVCF files")
        return [str(f) for f in gvcf_files]

    def read_chromosomes(self, reference_fai: str) -> List[str]:
        """
        读取染色体列表|Read chromosome list

        Args:
            reference_fai: 参考基因组索引文件|Reference genome index file

        Returns:
            染色体列表|List of chromosomes
        """
        self.logger.info(f"正在读取染色体信息...|Reading chromosome information...")

        try:
            with open(reference_fai, 'r') as f:
                chromosomes = [line.split('\t')[0] for line in f]

            self.logger.info(f"参考基因组共有 {len(chromosomes)} 条染色体/scaffold|"
                           f"Reference has {len(chromosomes)} chromosomes/scaffolds")
            return chromosomes
        except Exception as e:
            self.logger.error(f"读取染色体信息失败|Failed to read chromosome information: {e}")
            return []

    def check_faketime(self) -> bool:
        """
        检查faketime是否可用|Check if faketime is available

        Returns:
            faketime是否可用|Whether faketime is available
        """
        return shutil.which('faketime') is not None

    def run_command(self, command: str, description: str = None) -> bool:
        """
        执行命令并记录日志|Execute command and log output

        Args:
            command: 要执行的命令|Command to execute
            description: 命令描述|Command description

        Returns:
            是否成功|Whether successful
        """
        if description:
            self.logger.info(f"{description}|{description}")

        try:
            result = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                check=True
            )
            if result.stdout:
                self.logger.debug(f"命令输出|Command output: {result.stdout.strip()}")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command execution failed: {e}")
            if e.stderr:
                self.logger.error(f"错误信息|Error output: {e.stderr.strip()}")
            return False
        except Exception as e:
            self.logger.error(f"执行出错|Execution error: {e}")
            return False

    def ensure_reference_index(self, reference: str) -> bool:
        """
        确保参考基因组索引存在，不存在则自动创建|Ensure reference genome index exists, create if not

        Args:
            reference: 参考基因组文件路径|Reference genome file path

        Returns:
            是否成功|Whether successful
        """
        fai_file = f"{reference}.fai"

        if os.path.exists(fai_file):
            self.logger.info(f"参考基因组索引已存在|Reference index exists: {fai_file}")
            return True

        self.logger.warning(f"参考基因组索引不存在|Reference index not found: {fai_file}")
        self.logger.info(f"正在创建索引...|Creating index...")

        command = f"samtools faidx {reference}"
        if self.run_command(command, f"运行 samtools faidx|Running samtools faidx"):
            self.logger.info(f"索引创建成功|Index created successfully: {fai_file}")
            return True
        else:
            self.logger.error(f"索引创建失败|Failed to create index")
            return False

    def check_and_rebuild_vcf_index(self, vcf_file: str) -> bool:
        """
        检查VCF索引，如果索引过时则重建|Check VCF index, rebuild if outdated

        Args:
            vcf_file: VCF文件路径|VCF file path

        Returns:
            是否成功|Whether successful
        """
        tbi_file = f"{vcf_file}.tbi"

        # 检查索引是否存在
        if not os.path.exists(tbi_file):
            self.logger.warning(f"VCF索引不存在|VCF index not found: {tbi_file}")
            return self._rebuild_vcf_index(vcf_file)

        # 检查时间戳
        vcf_mtime = os.path.getmtime(vcf_file)
        tbi_mtime = os.path.getmtime(tbi_file)

        if tbi_mtime < vcf_mtime:
            self.logger.warning(
                f"VCF索引过时|VCF index outdated: {tbi_file} "
                f"(索引时间|Index time: {tbi_mtime} < VCF时间|VCF time: {vcf_mtime})"
            )
            self.logger.info(f"删除旧索引|Removing old index: {tbi_file}")
            try:
                os.remove(tbi_file)
            except Exception as e:
                self.logger.error(f"删除索引失败|Failed to remove index: {e}")
                return False
            return self._rebuild_vcf_index(vcf_file)

        return True

    def _rebuild_vcf_index(self, vcf_file: str) -> bool:
        """
        重建VCF索引|Rebuild VCF index

        Args:
            vcf_file: VCF文件路径|VCF file path

        Returns:
            是否成功|Whether successful
        """
        command = f"tabix -p vcf {vcf_file}"
        if self.run_command(command, f"重建VCF索引|Rebuilding VCF index"):
            self.logger.info(f"索引重建成功|Index rebuilt successfully: {vcf_file}.tbi")
            return True
        else:
            self.logger.error(f"索引重建失败|Failed to rebuild index")
            return False

    def scan_and_validate_gvcf_files(self, gvcf_dir: str, rebuild_indexes: bool = True) -> List[str]:
        """
        扫描并验证GVCF文件|Scan and validate GVCF files

        Args:
            gvcf_dir: GVCF目录|GVCF directory
            rebuild_indexes: 是否重建过时的索引|Whether to rebuild outdated indexes

        Returns:
            GVCF文件列表|List of GVCF files
        """
        self.logger.info(f"正在扫描 GVCF 文件...|Scanning GVCF files...")
        gvcf_files = sorted(Path(gvcf_dir).glob("*.g.vcf.gz"))

        if not gvcf_files:
            self.logger.error(f"未找到任何 *.g.vcf.gz 文件|No *.g.vcf.gz files found")
            return []

        valid_files = []
        for vcf_file in gvcf_files:
            vcf_str = str(vcf_file)

            # 检查并重建索引
            if rebuild_indexes:
                if not self.check_and_rebuild_vcf_index(vcf_str):
                    self.logger.warning(f"跳过文件|Skipping file: {vcf_str}")
                    continue

            valid_files.append(vcf_str)

        self.logger.info(f"找到 {len(valid_files)} 个有效 GVCF 文件|Found {len(valid_files)} valid GVCF files")
        return valid_files


class GTXJointWriter:
    """GTX Joint Calling命令写入器|GTX Joint Calling Command Writer"""

    def __init__(self, logger):
        """
        初始化写入器|Initialize writer

        Args:
            logger: 日志对象|Logger object
        """
        self.logger = logger

    def backup_old_script(self, script_path: str) -> bool:
        """
        备份旧脚本|Backup old script

        Args:
            script_path: 脚本路径|Script path

        Returns:
            是否成功备份|Whether backup was successful
        """
        if os.path.exists(script_path):
            backup = f"{script_path}.bak.{self._get_timestamp()}"
            try:
                shutil.move(script_path, backup)
                self.logger.warning(f"已备份旧文件为|Old script backed up as: {backup}")
                return True
            except Exception as e:
                self.logger.error(f"备份失败|Backup failed: {e}")
                return False
        return True

    def write_command(self, script_path: str, command: str) -> bool:
        """
        写入命令到脚本|Write command to script

        Args:
            script_path: 脚本路径|Script path
            command: 命令字符串|Command string

        Returns:
            是否成功写入|Whether write was successful
        """
        try:
            with open(script_path, 'a') as f:
                f.write(command + '\n')
            return True
        except Exception as e:
            self.logger.error(f"写入命令失败|Failed to write command: {e}")
            return False

    def make_executable(self, script_path: str) -> bool:
        """
        设置脚本为可执行|Make script executable

        Args:
            script_path: 脚本路径|Script path

        Returns:
            是否成功设置|Whether setting was successful
        """
        try:
            os.chmod(script_path, 0o755)
            return True
        except Exception as e:
            self.logger.error(f"设置可执行权限失败|Failed to set executable permission: {e}")
            return False

    def _get_timestamp(self) -> str:
        """获取时间戳|Get timestamp"""
        from datetime import datetime
        return datetime.now().strftime("%Y%m%d_%H%M%S")
