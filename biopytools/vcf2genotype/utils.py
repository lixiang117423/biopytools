"""
VCF基因型提取工具函数模块|VCF Genotype Extraction Utility Functions Module
"""

import gzip
import logging
import sys
import os
import time
from pathlib import Path
from typing import Optional


class VCFLogger:
    """VCF提取日志管理器|VCF Extraction Logger Manager"""

    def __init__(self, log_dir: str, verbose: bool = False, quiet: bool = False):
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.verbose = verbose
        self.quiet = quiet

        # 创建日志文件|Create log file
        timestamp = time.strftime('%Y%m%d_%H%M%S')
        logs_dir = self.log_dir / "99_logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        self.log_file = logs_dir / f"vcf_extraction_{timestamp}.log"

        # 配置logger|Configure logger
        self.logger = logging.getLogger(f"vcf_extraction_{timestamp}")
        self.logger.propagate = False  # 避免重复输出到root logger|Avoid duplicate output to root logger

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

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
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


class FileUtils:
    """文件工具类|File Utility Class"""

    @staticmethod
    def is_gzipped(file_path: str) -> bool:
        """检查文件是否为gzip压缩|Check if file is gzip compressed"""
        try:
            with gzip.open(file_path, 'rt') as f:
                f.readline()
            return True
        except Exception:
            return False

    @staticmethod
    def open_file(file_path: str, mode: str = 'r'):
        """智能打开文件（自动检测是否压缩）|Smart file opening (auto-detect compression)"""
        if FileUtils.is_gzipped(file_path):
            return gzip.open(file_path, mode + 't', encoding='utf-8')
        else:
            return open(file_path, mode, encoding='utf-8')


class GenotypeUtils:
    """基因型工具类|Genotype Utility Class"""

    @staticmethod
    def parse_genotype(gt_field: str) -> Optional[str]:
        """解析基因型字段|Parse genotype field"""
        if gt_field in ['.', './.', '.|.']:
            return None

        # 移除相位信息，只保留基因型|Remove phasing info, keep only genotype
        gt = gt_field.split(':')[0] if ':' in gt_field else gt_field

        return gt

    @staticmethod
    def is_biallelic(ref: str, alt: str) -> bool:
        """检查是否为双等位变异|Check if variant is biallelic"""
        # 如果ALT字段不包含逗号，则为双等位|If ALT field doesn't contain comma, it's biallelic
        return ',' not in alt

    @staticmethod
    def calculate_genotype_stats(genotypes: list) -> tuple:
        """计算基因型统计信息|Calculate genotype statistics"""
        valid_gts = [gt for gt in genotypes if gt is not None]

        if not valid_gts:
            return 0.0, 0.0

        homozygous_count = 0
        heterozygous_count = 0

        for gt in valid_gts:
            # 处理不同的分隔符|Handle different separators
            if '/' in gt:
                alleles = gt.split('/')
            elif '|' in gt:
                alleles = gt.split('|')
            else:
                continue

            if len(alleles) == 2:
                if alleles[0] == alleles[1]:
                    homozygous_count += 1
                else:
                    heterozygous_count += 1

        total = len(valid_gts)
        homozygous_ratio = homozygous_count / total if total > 0 else 0.0
        heterozygous_ratio = heterozygous_count / total if total > 0 else 0.0

        return homozygous_ratio, heterozygous_ratio


def check_dependencies():
    """检查依赖|Check dependencies"""
    optional_deps = {}

    try:
        import cyvcf2
        optional_deps['cyvcf2'] = True
    except ImportError:
        optional_deps['cyvcf2'] = False

    try:
        import pandas
        optional_deps['pandas'] = True
    except ImportError:
        optional_deps['pandas'] = False

    return optional_deps
