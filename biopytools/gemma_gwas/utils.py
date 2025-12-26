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
    required_commands = ['plink', 'awk']
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


def fix_fam_file(pheno_file: str, fam_file: str, logger: logging.Logger) -> bool:
    """
    修复FAM文件的表型列

    Args:
        pheno_file: 表型文件路径
        fam_file: FAM文件路径
        logger: 日志对象

    Returns:
        bool: 是否成功
    """
    try:
        # 提取第一个表型列
        pheno_values = []
        with open(pheno_file, 'r') as f:
            lines = f.readlines()
            if len(lines) < 2:
                logger.error("Phenotype file has no data rows")
                return False

            for line in lines[1:]:  # 跳过表头
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    pheno_values.append(fields[1])

        # 读取并修改FAM文件
        with open(fam_file, 'r') as f:
            fam_lines = f.readlines()

        # 修改第6列（表型列）
        with open(f"{fam_file}.tmp", 'w') as f:
            for i, line in enumerate(fam_lines):
                fields = line.strip().split()
                if len(fields) >= 6 and i < len(pheno_values):
                    fields[5] = pheno_values[i]
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

        # 跳过表头，提取第3列及以后的主成分
        with open(output_file, 'w') as f:
            for line in lines[1:]:
                fields = line.strip().split()
                # 从第3列开始（前两列是FID和IID）
                pc_values = '\t'.join(fields[2:])
                f.write(pc_values + '\n')

        logger.info(f"Covariate file prepared: {output_file}")
        return True

    except Exception as e:
        logger.error(f"Failed to prepare covariate file: {str(e)}")
        return False


def match_samples(pca_file: str, genotype_prefix: str, pheno_file: str,
                  output_prefix: str, threads: int, logger: logging.Logger) -> bool:
    """
    同步基因型、表型和PCA的样本

    Args:
        pca_file: PCA文件路径
        genotype_prefix: 基因型文件前缀
        pheno_file: 表型文件路径
        output_prefix: 输出文件前缀
        threads: 线程数
        logger: 日志对象

    Returns:
        bool: 是否成功
    """
    try:
        # 从PCA文件提取样本ID
        keep_file = f"{output_prefix}_keep_samples.txt"

        with open(pca_file, 'r') as f:
            lines = f.readlines()

        with open(keep_file, 'w') as f:
            for line in lines[1:]:  # 跳过表头
                fields = line.strip().split()
                if len(fields) >= 2:
                    f.write(f"{fields[0]}\t{fields[1]}\n")

        # 使用PLINK过滤样本
        cmd = f"plink --bfile {genotype_prefix} --keep {keep_file} " \
              f"--make-bed --out {output_prefix} " \
              f"--allow-extra-chr --threads {threads}"

        log_file = f"{output_prefix}_match.log"
        success, error = run_command(cmd, log_file, logger)

        if not success:
            logger.error(f"Failed to match samples: {error}")
            return False

        # 过滤表型文件
        keep_samples = set()
        with open(keep_file, 'r') as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 2:
                    keep_samples.add(fields[1])  # IID

        pheno_matched_file = f"{output_prefix}_pheno_matched.txt"
        with open(pheno_file, 'r') as f_in:
            with open(pheno_matched_file, 'w') as f_out:
                # 保留表头
                header = f_in.readline()
                f_out.write(header)

                # 过滤样本
                for line in f_in:
                    fields = line.strip().split('\t')
                    if fields[0] in keep_samples:
                        f_out.write(line)

        # 替换原表型文件
        os.replace(pheno_matched_file, pheno_file)

        # 清理临时文件
        os.remove(keep_file)

        logger.info("Samples matched successfully")
        return True

    except Exception as e:
        logger.error(f"Failed to match samples: {str(e)}")
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
