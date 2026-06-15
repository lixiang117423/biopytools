"""
覆盖度过滤核心计算模块|Coverage Filter Core Calculation Module
"""

import os
import subprocess
from pathlib import Path
from .utils import run_command


class CoverageCalculator:
    """覆盖度计算器|Coverage Calculator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def calculate_coverage(self):
        """计算覆盖度|Calculate coverage"""
        self.logger.info("开始计算覆盖度|Starting coverage calculation")

        coverage_file = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_coverage.txt")

        cmd = f"samtools coverage {self.config.bam_file} > {coverage_file}"

        if not run_command(cmd, self.logger, "覆盖度计算|Coverage calculation"):
            return None

        self.logger.info(f"覆盖度文件已生成|Coverage file generated: {coverage_file}")
        return coverage_file

    def classify_sequences(self, coverage_file):
        """分类序列质量|Classify sequence quality"""
        self.logger.info("开始分类序列质量|Starting sequence quality classification")

        # 高质量: 覆盖度≥high_cov|High quality: coverage≥high_cov
        high_quality_list = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_high_quality.list")
        cmd = (f"awk 'NR>1 && $6>={self.config.high_coverage_min} {{print $1}}' "
               f"{coverage_file} > {high_quality_list}")

        if not run_command(cmd, self.logger, "提取高质量序列|Extracting high quality sequences"):
            return None, None, None

        high_count = self._count_lines(high_quality_list)
        self.logger.info(f"高质量序列|High quality sequences: {high_count}")

        # 中等质量: 覆盖度medium_cov_min到high_cov|Medium quality: coverage medium_cov_min to high_cov
        medium_quality_list = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_medium_quality.list")
        cmd = (f"awk 'NR>1 && $6>={self.config.medium_coverage_min} "
               f"&& $6<{self.config.medium_coverage_max} "
               f"{{print $1}}' {coverage_file} > {medium_quality_list}")

        if not run_command(cmd, self.logger, "提取中等质量序列|Extracting medium quality sequences"):
            return None, None, None

        medium_count = self._count_lines(medium_quality_list)
        self.logger.info(f"中等质量序列|Medium quality sequences: {medium_count}")

        # 低质量: 覆盖度<medium_cov_min|Low quality: coverage<medium_cov_min
        low_quality_list = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_low_quality.list")
        cmd = (f"awk 'NR>1 && $6<{self.config.low_coverage_max} {{print $1}}' "
               f"{coverage_file} > {low_quality_list}")

        if not run_command(cmd, self.logger, "提取低质量序列|Extracting low quality sequences"):
            return None, None, None

        low_count = self._count_lines(low_quality_list)
        self.logger.info(f"低质量序列|Low quality sequences: {low_count}")

        return high_quality_list, medium_quality_list, low_quality_list

    def extract_sequences(self, high_list, medium_list, low_list):
        """提取序列|Extract sequences"""
        self.logger.info("开始提取序列|Starting sequence extraction")

        # 提取高质量序列|Extract high quality sequences
        high_fasta = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_high_quality.fa")
        cmd = f"seqtk subseq {self.config.fasta_file} {high_list} > {high_fasta}"

        if not run_command(cmd, self.logger, "提取高质量序列|Extracting high quality sequences"):
            return False

        # 提取中等质量序列|Extract medium quality sequences
        medium_fasta = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_medium_quality.fa")
        cmd = f"seqtk subseq {self.config.fasta_file} {medium_list} > {medium_fasta}"

        if not run_command(cmd, self.logger, "提取中等质量序列|Extracting medium quality sequences"):
            return False

        # 提取低质量序列|Extract low quality sequences
        low_fasta = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_low_quality.fa")
        cmd = f"seqtk subseq {self.config.fasta_file} {low_list} > {low_fasta}"

        if not run_command(cmd, self.logger, "提取低质量序列|Extracting low quality sequences"):
            return False

        self.logger.info("序列提取完成|Sequence extraction completed")
        return True

    def generate_statistics(self):
        """生成统计信息|Generate statistics"""
        self.logger.info("开始生成统计信息|Starting statistics generation")

        stats = {}

        # 原始序列统计|Original sequence statistics
        cmd = f"seqkit stats -j {self.config.threads} {self.config.fasta_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            stats['original'] = result.stdout
            self.logger.info("原始序列统计完成|Original sequence statistics completed")
        else:
            self.logger.warning("原始序列统计失败|Original sequence statistics failed")
            stats['original'] = "统计失败|Statistics failed"

        # 高质量序列统计|High quality sequence statistics
        high_fasta = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_high_quality.fa")
        if os.path.exists(high_fasta):
            cmd = f"seqkit stats -j {self.config.threads} {high_fasta}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                stats['high_quality'] = result.stdout
                self.logger.info("高质量序列统计完成|High quality sequence statistics completed")
            else:
                stats['high_quality'] = "统计失败|Statistics failed"
        else:
            stats['high_quality'] = "文件不存在|File not found"

        # 中等质量序列统计|Medium quality sequence statistics
        medium_fasta = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_medium_quality.fa")
        if os.path.exists(medium_fasta):
            cmd = f"seqkit stats -j {self.config.threads} {medium_fasta}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                stats['medium_quality'] = result.stdout
                self.logger.info("中等质量序列统计完成|Medium quality sequence statistics completed")
            else:
                stats['medium_quality'] = "统计失败|Statistics failed"
        else:
            stats['medium_quality'] = "文件不存在|File not found"

        # 低质量序列统计|Low quality sequence statistics
        low_fasta = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_low_quality.fa")
        if os.path.exists(low_fasta):
            cmd = f"seqkit stats -j {self.config.threads} {low_fasta}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                stats['low_quality'] = result.stdout
                self.logger.info("低质量序列统计完成|Low quality sequence statistics completed")
            else:
                stats['low_quality'] = "统计失败|Statistics failed"
        else:
            stats['low_quality'] = "文件不存在|File not found"

        return stats

    def _count_lines(self, file_path):
        """统计文件行数|Count file lines"""
        try:
            with open(file_path, 'r') as f:
                return sum(1 for _ in f)
        except FileNotFoundError:
            return 0
