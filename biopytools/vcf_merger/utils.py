"""
VCF合并工具实用函数模块|VCF Merger Utility Functions Module
"""

import logging
import sys
import subprocess
import re
from pathlib import Path
from typing import List, Optional, Dict
from collections import defaultdict


class VCFMergerLogger:
    """VCF合并工具日志管理器|VCF Merger Logger Manager"""

    def __init__(self, log_file: Optional[str] = None, level: str = "INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path (optional)
            level: 日志级别|Log level (default: INFO)
        """
        self.logger = logging.getLogger("VCFMerger")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()

        # 日志格式|Log format
        formatter = logging.Formatter(
            "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )

        # stdout handler - INFO级别及以下|stdout handler - INFO and below
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上级别|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler(如果指定)|File handler (if specified)
        if log_file:
            file_handler = logging.FileHandler(log_file, encoding="utf-8")
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

        # 设置日志级别|Set log level
        level_map = {
            "DEBUG": logging.DEBUG,
            "INFO": logging.INFO,
            "WARNING": logging.WARNING,
            "ERROR": logging.ERROR,
            "CRITICAL": logging.CRITICAL
        }
        self.logger.setLevel(level_map.get(level.upper(), logging.INFO))

    def get_logger(self):
        """获取logger对象|Get logger object"""
        return self.logger


class VCFFileFinder:
    """VCF文件查找器|VCF File Finder"""

    @staticmethod
    def find_vcf_files(input_dir: Path, pattern: str, logger: logging.Logger) -> List[Path]:
        """
        查找所有VCF文件|Find all VCF files

        Args:
            input_dir: 输入目录|Input directory
            pattern: 文件匹配模式|File pattern
            logger: 日志对象|Logger object

        Returns:
            VCF文件列表|List of VCF files

        Raises:
            FileNotFoundError: 未找到VCF文件时|When no VCF files found
        """
        if not input_dir.exists():
            raise FileNotFoundError(
                f"输入目录不存在|Input directory does not exist: {input_dir}"
            )

        vcf_files = list(input_dir.glob(pattern))

        if not vcf_files:
            raise FileNotFoundError(
                f"未找到匹配的VCF文件|No VCF files found matching pattern: {pattern}"
            )

        logger.info(f"找到 {len(vcf_files)} 个VCF文件|Found {len(vcf_files)} VCF files")
        return vcf_files


class ChromosomeExtractor:
    """染色体编号提取器|Chromosome Number Extractor"""

    @staticmethod
    def extract_chromosome(filename: str, logger: logging.Logger) -> Optional[str]:
        """
        从文件名中提取染色体编号|Extract chromosome number from filename

        支持格式|Supported formats:
        - Chr19_1-20000000.joint.vcf.gz -> Chr19
        - chr19_1-20000000.joint.vcf.gz -> chr19
        - 19_1-20000000.joint.vcf.gz -> 19

        Args:
            filename: 文件名|Filename
            logger: 日志对象|Logger object

        Returns:
            染色体编号|Chromosome ID or None if not found
        """
        # 匹配Chr/chr开头或纯数字的染色体编号（只匹配到第一个下划线之前）
        # Match chromosome numbers starting with Chr/chr or pure numbers (stop before first underscore)
        match = re.match(r"^([Cc]hr)?(\w+)(?=(_|$))", filename)

        if match:
            prefix = match.group(1) or ""
            chr_num = match.group(2)
            chr_id = f"{prefix}{chr_num}"
            logger.debug(
                f"从文件名提取染色体|Extracted chromosome from filename: "
                f"{filename} -> {chr_id}"
            )
            return chr_id

        logger.warning(
            f"无法从文件名提取染色体编号|Cannot extract chromosome from filename: {filename}"
        )
        return None


class VCFGroupByChromosome:
    """按染色体分组VCF文件|Group VCF files by chromosome"""

    @staticmethod
    def group_vcf_by_chromosome(
        vcf_files: List[Path],
        logger: logging.Logger
    ) -> Dict[str, List[Path]]:
        """
        按染色体分组VCF文件|Group VCF files by chromosome

        Args:
            vcf_files: VCF文件列表|List of VCF files
            logger: 日志对象|Logger object

        Returns:
            按染色体分组的VCF文件字典|Dictionary of VCF files grouped by chromosome
        """
        chr_groups = defaultdict(list)

        for vcf_file in vcf_files:
            chr_id = ChromosomeExtractor.extract_chromosome(vcf_file.name, logger)

            if chr_id:
                chr_groups[chr_id].append(vcf_file)
                logger.debug(
                    f"文件 {vcf_file.name} 分配到染色体 {chr_id}|"
                    f"File {vcf_file.name} assigned to chromosome {chr_id}"
                )

        return chr_groups


class BCFToolsChecker:
    """bcftools可用性检查器|BCFTools Availability Checker"""

    @staticmethod
    def check_bcftools(logger: logging.Logger) -> bool:
        """
        检查bcftools是否可用|Check if bcftools is available

        Args:
            logger: 日志对象|Logger object

        Returns:
            bcftools是否可用|Whether bcftools is available
        """
        try:
            result = subprocess.run(
                ["bcftools", "--version"],
                capture_output=True,
                text=True,
                check=False
            )

            if result.returncode == 0:
                version_line = result.stdout.split("\n")[0]
                logger.info(f"bcftools已安装|bcftools is available: {version_line}")
                return True

        except FileNotFoundError:
            pass

        logger.error(
            "未找到bcftools，请先安装|bcftools not found, please install first\n"
            "安装命令|Install command: conda install -c bioconda bcftools"
        )
        return False


class VCFMerger:
    """VCF文件合并器|VCF File Merger using bcftools"""

    @staticmethod
    def merge_vcf_files(
        vcf_list: List[Path],
        output_file: Path,
        threads: int,
        logger: logging.Logger
    ) -> bool:
        """
        使用bcftools concat合并VCF文件|Merge VCF files using bcftools concat

        Args:
            vcf_list: VCF文件列表|List of VCF files
            output_file: 输出文件路径|Output file path
            threads: 线程数|Number of threads
            logger: 日志对象|Logger object

        Returns:
            是否成功合并|Whether merge was successful
        """
        logger.info(
            f"合并 {len(vcf_list)} 个文件到 {output_file.name}|"
            f"Merging {len(vcf_list)} files to {output_file.name}"
        )

        # 构建输入文件列表|Build input file list
        input_list_file = output_file.parent / f"{output_file.stem}_filelist.txt"

        with open(input_list_file, "w") as f:
            for vcf in sorted(vcf_list):
                f.write(f"{vcf}\n")

        try:
            # bcftools concat命令|bcftools concat command
            cmd = [
                "bcftools",
                "concat",
                "-f", str(input_list_file),
                "-O", "z",  # 输出压缩的VCF|Output compressed VCF
                "-o", str(output_file),
                "--threads", str(threads)
            ]

            logger.debug(f"执行命令|Executing command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )

            if result.returncode != 0:
                logger.error(f"合并失败|Merge failed: {result.stderr}")
                return False

            logger.info(f"成功创建 {output_file.name}|Successfully created {output_file.name}")

            # 删除临时文件列表|Delete temporary file list
            input_list_file.unlink()

            return True

        except Exception as e:
            logger.error(f"合并过程中出错|Error during merge: {e}")
            return False


class VCFIndexer:
    """VCF索引创建器|VCF Index Creator using bcftools"""

    @staticmethod
    def index_vcf_file(vcf_file: Path, logger: logging.Logger) -> bool:
        """
        使用bcftools index创建索引|Create index using bcftools index

        Args:
            vcf_file: VCF文件路径|VCF file path
            logger: 日志对象|Logger object

        Returns:
            是否成功创建索引|Whether index creation was successful
        """
        logger.info(f"正在创建索引|Creating index for: {vcf_file.name}")

        try:
            cmd = ["bcftools", "index", str(vcf_file)]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )

            if result.returncode != 0:
                logger.warning(f"索引创建失败|Index creation failed: {result.stderr}")
                return False

            logger.info(f"索引创建成功|Index created successfully")
            return True

        except Exception as e:
            logger.warning(f"索引创建过程中出错|Error during indexing: {e}")
            return False
