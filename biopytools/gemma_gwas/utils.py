#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GEMMA GWAS Utility Functions
工具函数模块
Author: BioTools Development Team
Version: 1.0.0
"""

import os
import subprocess
import logging
from typing import Tuple, List
from pathlib import Path


def run_command(cmd: str, log_file: str, logger: logging.Logger) -> Tuple[bool, str]:
    """
    运行shell命令并记录日志

    Args:
        cmd: 要执行的命令
        log_file: 日志文件路径
        logger: 日志对象

    Returns:
        tuple: (是否成功, 错误消息)
    """
    try:
        with open(log_file, 'w') as f:
            result = subprocess.run(
                cmd,
                shell=True,
                stdout=f,
                stderr=subprocess.STDOUT,
                text=True
            )

        if result.returncode != 0:
            return False, f"Command failed with return code {result.returncode}"

        return True, ""

    except Exception as e:
        return False, f"Exception: {str(e)}"


def check_dependencies(logger: logging.Logger) -> Tuple[bool, str]:
    """
    检查依赖程序是否安装

    Args:
        logger: 日志对象

    Returns:
        tuple: (是否都安装, 缺失的程序列表)
    """
    required_commands = ['plink', 'bcftools', 'awk']
    missing = []

    for cmd in required_commands:
        try:
            result = subprocess.run(
                ['which', cmd],
                capture_output=True,
                text=True
            )
            if result.returncode != 0:
                missing.append(cmd)
        except Exception:
            missing.append(cmd)

    if missing:
        return False, ", ".join(missing)

    return True, ""


def extract_vcf_samples(vcf_file: str, logger: logging.Logger = None) -> set:
    """
    从VCF文件中提取样本ID（不需要完全转换）

    Args:
        vcf_file: VCF文件路径
        logger: 日志对象

    Returns:
        set: 样本ID集合
    """
    try:
        # 使用bcftools查询样本ID
        cmd = f"bcftools query -l {vcf_file}"
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            if logger:
                logger.error(f"Failed to extract samples from VCF: {result.stderr}")
            return set()

        samples = set()
        for line in result.stdout.strip().split('\n'):
            if line:
                samples.add(line)

        return samples

    except Exception as e:
        if logger:
            logger.error(f"Failed to extract VCF samples: {str(e)}")
        return set()


def count_lines(filepath: str) -> int:
    """
    统计文件行数

    Args:
        filepath: 文件路径

    Returns:
        int: 行数
    """
    try:
        with open(filepath, 'r') as f:
            return sum(1 for _ in f)
    except Exception:
        return 0


def count_columns(filepath: str, delimiter: str = '\t') -> int:
    """
    统计文件列数（基于第一行）

    Args:
        filepath: 文件路径
        delimiter: 分隔符

    Returns:
        int: 列数
    """
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline().strip()
            return len(first_line.split(delimiter))
    except Exception:
        return 0


def fix_fam_file(pheno_file: str, fam_file: str, logger: logging.Logger,
                 has_header: bool = False) -> bool:
    """
    修复FAM文件的表型列

    Args:
        pheno_file: 表型文件路径
        fam_file: FAM文件路径
        logger: 日志对象
        has_header: 表型文件是否有表头

    Returns:
        bool: 是否成功
    """
    try:
        # 提取第一个表型列
        pheno_values = []
        pheno_dict = {}  # 样本ID到表型值的映射

        with open(pheno_file, 'r') as f:
            lines = f.readlines()
            if len(lines) < 1:
                logger.error("Phenotype file is empty")
                return False

            start_idx = 1 if has_header else 0

            for line in lines[start_idx:]:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    sample_id = fields[0]
                    pheno_value = fields[1]
                    pheno_dict[sample_id] = pheno_value

        # 读取并修改FAM文件
        with open(fam_file, 'r') as f:
            fam_lines = f.readlines()

        # 修改第6列（表型列）
        with open(f"{fam_file}.tmp", 'w') as f:
            for line in fam_lines:
                fields = line.strip().split()
                if len(fields) >= 2:
                    sample_id = fields[1]  # IID
                    if sample_id in pheno_dict:
                        fields[5] = pheno_dict[sample_id]
                f.write('\t'.join(fields) + '\n')

        # 替换原文件
        os.replace(f"{fam_file}.tmp", fam_file)

        logger.info("FAM file phenotype column updated")
        return True

    except Exception as e:
        logger.error(f"Failed to fix FAM file: {str(e)}")
        return False


def prepare_phenotype_files(pheno_file: str, output_dir: str,
                            logger: logging.Logger) -> Tuple[bool, List[str]]:
    """
    准备表型文件（带表头和不带表头）

    Args:
        pheno_file: 表型文件路径
        output_dir: 输出目录
        logger: 日志对象

    Returns:
        tuple: (是否成功, [pheno_with_header, pheno_no_header, pheno_names])
    """
    try:
        base_name = Path(pheno_file).stem

        # 读取表头
        with open(pheno_file, 'r') as f:
            header = f.readline().strip()
            pheno_names = header.split('\t')

        pheno_with_header = os.path.join(output_dir, f"{base_name}_with_header.txt")
        pheno_no_header = os.path.join(output_dir, f"{base_name}_no_header.txt")

        # 带表头的文件（包含所有行）
        with open(pheno_file, 'r') as f:
            with open(pheno_with_header, 'w') as out:
                out.write(f.read())

        # 不带表头的文件（跳过第一行）
        with open(pheno_file, 'r') as f:
            next(f)  # 跳过表头
            with open(pheno_no_header, 'w') as out:
                out.write(f.read())

        logger.info(f"Prepared phenotype files:")
        logger.info(f"  With header: {pheno_with_header}")
        logger.info(f"  No header: {pheno_no_header}")
        logger.info(f"  Phenotypes: {', '.join(pheno_names)}")

        return True, [pheno_with_header, pheno_no_header, pheno_names]

    except Exception as e:
        logger.error(f"Failed to prepare phenotype files: {str(e)}")
        return False, []


def find_common_samples_and_filter(
    vcf_file: str,
    pheno_file_no_header: str,
    pheno_with_header: str,
    output_prefix: str,
    threads: int,
    logger: logging.Logger
) -> Tuple[bool, str, str]:
    """
    找到VCF和表型的样本交集，并使用bcftools从VCF筛选样本

    Args:
        vcf_file: 原始VCF文件路径
        pheno_file_no_header: 表型文件路径（无表头）
        pheno_with_header: 表型文件路径（带表头）
        output_prefix: 输出文件前缀
        threads: 线程数
        logger: 日志对象

    Returns:
        tuple: (是否成功, 过滤后的表型文件路径(无表头), 样本数量)
    """
    try:
        # 1. 从VCF文件读取样本ID（使用bcftools，快速）
        vcf_samples = extract_vcf_samples(vcf_file, logger)

        if len(vcf_samples) == 0:
            logger.error("No samples found in VCF file!")
            return False, "", 0

        logger.info(f"VCF samples: {len(vcf_samples)}")

        # 2. 从表型文件读取样本ID
        pheno_samples = set()
        with open(pheno_file_no_header, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 1:
                    pheno_samples.add(fields[0])

        logger.info(f"Phenotype samples: {len(pheno_samples)}")

        # 3. 找到交集
        common_samples = vcf_samples.intersection(pheno_samples)
        logger.info(f"Common samples: {len(common_samples)}")

        if len(common_samples) == 0:
            logger.error("No common samples found between VCF and phenotype files!")
            return False, "", 0

        # 4. 使用bcftools从VCF筛选样本
        keep_file = f"{output_prefix}_keep.txt"
        vcf_filtered = f"{output_prefix}_filtered.vcf.gz"

        # 写入样本列表（每行一个样本ID）
        with open(keep_file, 'w') as f:
            for sample_id in sorted(common_samples):
                f.write(f"{sample_id}\n")

        logger.info(f"Filtering VCF from {len(vcf_samples)} to {len(common_samples)} samples...")

        # 使用bcftools筛选样本（多线程）
        cmd = f"bcftools view -S {keep_file} -Oz -o {vcf_filtered} {vcf_file} --threads {threads}"

        log_file = f"{output_prefix}_bcftools.log"
        success, error = run_command(cmd, log_file, logger)

        if not success:
            logger.error(f"Failed to filter VCF with bcftools: {error}")
            return False, "", 0

        logger.info("VCF filtering completed with bcftools")

        # 5. 将筛选后的VCF转换为PLINK格式
        logger.info("Converting filtered VCF to PLINK format...")
        cmd = f"plink --vcf {vcf_filtered} " \
              f"--make-bed --out {output_prefix} " \
              f"--allow-extra-chr --double-id --threads {threads}"

        log_file = f"{output_prefix}_plink_convert.log"
        success, error = run_command(cmd, log_file, logger)

        if not success:
            logger.error(f"Failed to convert filtered VCF to PLINK: {error}")
            return False, "", 0

        logger.info("PLINK conversion completed")

        # 6. 过滤表型文件（生成带表头和不带表头两个版本）
        pheno_filtered_with_header = f"{output_prefix}_pheno_filtered.txt"
        pheno_filtered_no_header = f"{output_prefix}_pheno_filtered_no_header.txt"

        with open(pheno_with_header, 'r') as f_in:
            with open(pheno_filtered_with_header, 'w') as f_out1, \
                 open(pheno_filtered_no_header, 'w') as f_out2:
                # 读取并写入表头（只写入带表头的文件）
                header = f_in.readline()
                f_out1.write(header)

                # 过滤样本
                for line in f_in:
                    fields = line.strip().split('\t')
                    if len(fields) >= 1 and fields[0] in common_samples:
                        f_out1.write(line)  # 写入带表头的文件
                        f_out2.write(line)  # 写入不带表头的文件

        # 7. 清理临时文件
        os.remove(keep_file)
        # 保留筛选后的VCF文件，可能有用

        logger.info(f"Filtered to {len(common_samples)} common samples")
        logger.info(f"Filtered phenotype files:")
        logger.info(f"  With header: {pheno_filtered_with_header}")
        logger.info(f"  No header: {pheno_filtered_no_header}")
        logger.info(f"Filtered VCF: {vcf_filtered}")

        return True, pheno_filtered_no_header, len(common_samples)

    except Exception as e:
        logger.error(f"Failed to find common samples: {str(e)}")
        return False, "", 0


def prepare_covariate_file(pca_file: str, output_file: str,
                          logger: logging.Logger) -> bool:
    """
    准备协变量文件（从PCA结果提取主成分）

    Args:
        pca_file: PCA结果文件
        output_file: 输出文件路径
        logger: 日志对象

    Returns:
        bool: 是否成功
    """
    try:
        with open(pca_file, 'r') as f:
            lines = f.readlines()

        # PCA eigenvec文件没有表头，不跳过第一行，提取第3列及以后的主成分
        with open(output_file, 'w') as f:
            for line in lines:
                fields = line.strip().split()
                # 从第3列开始（前两列是FID和IID）
                pc_values = '\t'.join(fields[2:])
                f.write(pc_values + '\n')

        logger.info(f"Covariate file prepared: {output_file} with {len(lines)} samples")
        return True

    except Exception as e:
        logger.error(f"Failed to prepare covariate file: {str(e)}")
        return False


def count_significant_snps(assoc_file: str, thresholds: List[float]) -> List[int]:
    """
    统计不同显著性阈值下的SNP数量

    Args:
        assoc_file: 关联分析结果文件
        thresholds: p值阈值列表

    Returns:
        list: 各阈值下的显著SNP数量
    """
    try:
        counts = [0] * len(thresholds)

        with open(assoc_file, 'r') as f:
            next(f)  # 跳过表头
            for line in f:
                fields = line.strip().split()
                if len(fields) < 1:
                    continue
                try:
                    p_value = float(fields[-1])  # 假设最后一列是p值
                    for i, threshold in enumerate(thresholds):
                        if p_value < threshold:
                            counts[i] += 1
                except ValueError:
                    continue

        return counts

    except Exception:
        return [0] * len(thresholds)
