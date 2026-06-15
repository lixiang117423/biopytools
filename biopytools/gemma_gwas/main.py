#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GEMMA GWAS Batch Analysis
GEMMA全基因组关联分析批量分析主流程
Author: BioTools Development Team
Version: 1.0.0
"""

import argparse
import logging
import sys
import os
import time
from pathlib import Path
from typing import List, Tuple

from .config import AnalysisConfig
from . import utils


# ==================== 日志配置 ====================

def setup_logger(name: str, log_file: str = None,
                level: int = logging.INFO) -> logging.Logger:
    """
    配置标准化的日志系统

    Args:
        name: logger名称
        log_file: 日志文件路径(可选)
        level: 日志级别

    Returns:
        logger对象
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

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
    logger.addHandler(stdout_handler)

    # stderr handler - WARNING及以上级别
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.WARNING)
    stderr_handler.setFormatter(formatter)
    logger.addHandler(stderr_handler)

    # 文件handler(如果指定)
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


# ==================== 分析步骤 ====================

def step1_prepare_phenotype(config: AnalysisConfig,
                            logger: logging.Logger) -> Tuple[bool, List[str]]:
    """
    步骤1: 准备表型文件

    Args:
        config: 分析配置
        logger: 日志对象

    Returns:
        tuple: (是否成功, [pheno_with_header, pheno_no_header, pheno_names])
    """
    logger.info("=" * 60)
    logger.info("STEP 1/6: Prepare Phenotype Files")
    logger.info("=" * 60)

    # 将表型文件路径转换为绝对路径
    pheno_abs = os.path.abspath(config.pheno) if not os.path.isabs(config.pheno) else config.pheno

    n_cols = utils.count_columns(pheno_abs)
    logger.info(f"Phenotype file columns: {n_cols} (including sample ID)")
    logger.info(f"Number of phenotypes: {n_cols - 1}")

    success, result = utils.prepare_phenotype_files(
        pheno_abs, config.outdir, logger
    )

    if not success:
        logger.error("Failed to prepare phenotype files")
        return False, []

    pheno_with_header, pheno_no_header, pheno_names = result

    # 转换为绝对路径以便在后续步骤中使用
    pheno_with_header_abs = os.path.abspath(pheno_with_header)
    pheno_no_header_abs = os.path.abspath(pheno_no_header)

    n_samples = utils.count_lines(pheno_no_header_abs)
    logger.info(f"Phenotype samples: {n_samples}")
    logger.info("Phenotype files prepared")
    logger.info("")

    return True, [pheno_with_header_abs, pheno_no_header_abs, pheno_names]


def step2_filter_common_samples(config: AnalysisConfig,
                                  pheno_with_header: str,
                                  pheno_no_header: str,
                                  logger: logging.Logger) -> Tuple[bool, str, int]:
    """
    步骤2: 从VCF和表型找到样本交集，并使用bcftools筛选

    Args:
        config: 分析配置
        pheno_with_header: 表型文件路径（带表头）
        pheno_no_header: 表型文件路径（无表头）
        logger: 日志对象

    Returns:
        tuple: (是否成功, 过滤后的表型文件路径, 样本数量)
    """
    logger.info("=" * 60)
    logger.info("STEP 2/6: Filter to Common Samples (using bcftools)")
    logger.info("=" * 60)

    # 切换到输出目录
    original_dir = os.getcwd()
    os.chdir(config.outdir)

    # 获取原始VCF文件的绝对路径
    vcf_abs = os.path.abspath(config.vcf) if not os.path.isabs(config.vcf) else config.vcf

    success, pheno_filtered, n_samples = utils.find_common_samples_and_filter(
        vcf_abs,
        pheno_no_header,
        pheno_with_header,
        "genotype",
        config.threads,
        logger
    )

    # 恢复原工作目录
    os.chdir(original_dir)

    if not success:
        logger.error("Failed to filter to common samples")
        return False, "", 0

    logger.info("Common samples filtering completed")
    logger.info("")

    # 返回绝对路径
    return True, os.path.abspath(pheno_filtered), n_samples


def step3_quality_control(config: AnalysisConfig,
                         logger: logging.Logger) -> bool:
    """
    步骤3: 质量控制（可选）

    Args:
        config: 分析配置
        logger: 日志对象

    Returns:
        bool: 是否成功
    """
    logger.info("=" * 60)
    logger.info("STEP 3/6: Quality Control")
    logger.info("=" * 60)

    # 切换到输出目录
    original_dir = os.getcwd()
    os.chdir(config.outdir)

    if not config.plink_qc.enable:
        logger.info("Quality control skipped")
        logger.info("")
        # 恢复原工作目录
        os.chdir(original_dir)
        return True

    # 构建质控命令
    qc_args = ' '.join(config.get_plink_qc_args())
    cmd = f"plink --bfile genotype " \
          f"--allow-extra-chr {qc_args} " \
          f"--make-bed --out genotype_qc " \
          f"--threads {config.threads}"

    log_file = "plink_qc.log"
    success, error = utils.run_command(cmd, log_file, logger)

    if not success:
        logger.error(f"Quality control failed: {error}")
        os.chdir(original_dir)
        return False

    # 统计信息
    bed_file = "genotype_qc.bed"
    if not os.path.exists(bed_file):
        logger.error("Output files not found")
        os.chdir(original_dir)
        return False

    n_samples_raw = utils.count_lines("genotype.fam")
    n_snps_raw = utils.count_lines("genotype.bim")
    n_samples_qc = utils.count_lines("genotype_qc.fam")
    n_snps_qc = utils.count_lines("genotype_qc.bim")

    logger.info(f"Quality control completed")
    logger.info(f"Samples after QC: {n_samples_qc} "
               f"(removed: {n_samples_raw - n_samples_qc})")
    logger.info(f"SNPs after QC: {n_snps_qc} "
               f"(removed: {n_snps_raw - n_snps_qc})")
    logger.info("")

    # 重命名genotype_qc为genotype
    for suffix in ['bed', 'bim', 'fam']:
        src = f"genotype_qc.{suffix}"
        dst = f"genotype.{suffix}"
        if os.path.exists(src):
            os.replace(src, dst)

    # 恢复原工作目录
    os.chdir(original_dir)

    return True


def step3b_fix_fam_file(config: AnalysisConfig,
                         pheno_filtered_no_header: str,
                         logger: logging.Logger) -> bool:
    """
    步骤3b: 修复FAM文件的表型列

    Args:
        config: 分析配置
        pheno_filtered_no_header: 过滤后的表型文件路径（无表头）
        logger: 日志对象

    Returns:
        bool: 是否成功
    """
    logger.info("Fixing FAM file phenotype column...")

    # 切换到输出目录
    original_dir = os.getcwd()
    os.chdir(config.outdir)

    fam_file = "genotype.fam"
    success = utils.fix_fam_file(pheno_filtered_no_header, fam_file, logger)

    # 恢复原工作目录
    os.chdir(original_dir)

    if success:
        logger.info("FAM file updated")
    else:
        logger.error("Failed to fix FAM file")

    logger.info("")
    return success


def step4_pca_analysis(config: AnalysisConfig,
                      pheno_no_header: str,
                      logger: logging.Logger) -> Tuple[bool, str]:
    """
    步骤4: PCA分析

    Args:
        config: 分析配置
        pheno_no_header: 表型文件路径（无表头）
        logger: 日志对象

    Returns:
        tuple: (是否成功, 协变量文件路径)
    """
    logger.info("=" * 60)
    logger.info(f"STEP 4/6: PCA Analysis (top {config.n_pca} components)")
    logger.info("=" * 60)

    # 切换到输出目录
    original_dir = os.getcwd()
    os.chdir(config.outdir)

    cmd = f"plink --bfile genotype " \
          f"--pca {config.n_pca} " \
          f"--out pca --allow-extra-chr " \
          f"--threads {config.threads}"

    log_file = "plink_pca.log"
    success, error = utils.run_command(cmd, log_file, logger)

    if not success:
        logger.error(f"PCA calculation failed: {error}")
        os.chdir(original_dir)
        return False, ""

    pca_file = "pca.eigenvec"
    if not os.path.exists(pca_file):
        logger.error("PCA output file not found")
        os.chdir(original_dir)
        return False, ""

    logger.info("PCA calculation completed")

    # 检查样本数
    pca_samples = utils.count_lines(pca_file)
    geno_samples = utils.count_lines("genotype.fam")

    logger.info(f"PCA samples: {pca_samples}")
    logger.info(f"Genotype samples: {geno_samples}")

    # 准备协变量文件
    covariate_file = "covariate.txt"
    success = utils.prepare_covariate_file(pca_file, covariate_file, logger)

    if not success:
        os.chdir(original_dir)
        return False, ""

    covariate_samples = utils.count_lines(covariate_file)
    logger.info(f"Covariate samples: {covariate_samples}")
    logger.info("Covariate file prepared")
    logger.info("")

    # 获取绝对路径返回
    covariate_file_abs = os.path.abspath(covariate_file)

    # 恢复原工作目录
    os.chdir(original_dir)

    return True, covariate_file_abs


def step5_kinship_matrix(config: AnalysisConfig,
                        pheno_no_header: str,
                        covariate_file: str,
                        logger: logging.Logger) -> bool:
    """
    步骤5: 计算亲缘关系矩阵

    Args:
        config: 分析配置
        pheno_no_header: 表型文件路径
        covariate_file: 协变量文件路径
        logger: 日志对象

    Returns:
        bool: 是否成功
    """
    logger.info("=" * 60)
    logger.info("STEP 5/6: Calculate Kinship Matrix")
    logger.info("=" * 60)

    # 最终样本数检查
    geno_file = os.path.join(config.outdir, "genotype.fam")
    final_geno = utils.count_lines(geno_file)
    final_pheno = utils.count_lines(pheno_no_header)
    final_covar = utils.count_lines(covariate_file)

    logger.info("Final sample count check:")
    logger.info(f"  Genotype: {final_geno}")
    logger.info(f"  Phenotype: {final_pheno}")
    logger.info(f"  Covariate: {final_covar}")

    if final_geno != final_pheno or final_geno != final_covar:
        logger.error("Sample count mismatch!")
        return False

    logger.info("Sample counts consistent")
    logger.info("")

    # 切换到输出目录
    original_dir = os.getcwd()
    os.chdir(config.outdir)

    # 构建GEMMA命令（GEMMA会自动创建output目录）
    qc_args = ' '.join(config.get_gemma_qc_args())
    cmd = f"{config.gemma_path} -bfile genotype " \
          f"-gk {config.gemma.gk_method} " \
          f"{qc_args} " \
          f"-o kinship"

    logger.info(f"Running command: {cmd}")

    log_file = "gemma_kinship.log"
    success, error = utils.run_command(cmd, log_file, logger)

    if not success:
        logger.error(f"Kinship matrix calculation failed: {error}")
        os.chdir(original_dir)
        return False

    # GEMMA会自动创建output目录
    kinship_file = "output/kinship.cXX.txt"
    if not os.path.exists(kinship_file):
        logger.error(f"Kinship matrix output file not found: {kinship_file}")
        os.chdir(original_dir)
        return False

    logger.info("Kinship matrix calculation completed")
    logger.info("")

    # 恢复原工作目录
    os.chdir(original_dir)

    return True


def step6_gwas_analysis(config: AnalysisConfig,
                       pheno_no_header: str,
                       covariate_file: str,
                       pheno_names: List[str],
                       logger: logging.Logger) -> bool:
    """
    步骤6: GWAS分析

    Args:
        config: 分析配置
        pheno_no_header: 表型文件路径
        covariate_file: 协变量文件路径
        pheno_names: 表型名称列表
        logger: 日志对象

    Returns:
        bool: 是否成功
    """
    logger.info("=" * 60)
    logger.info("STEP 6/6: GWAS Analysis")
    logger.info("=" * 60)
    logger.info("")

    # 切换到输出目录
    original_dir = os.getcwd()
    os.chdir(config.outdir)

    n_phenotypes = len(pheno_names) - 1  # 减去样本ID列
    gemma_output_dir = config.get_gemma_output_dir()

    # 分析每个表型
    for i in range(1, len(pheno_names)):  # 跳过第一列（样本ID）
        pheno_name = pheno_names[i]
        pheno_col = i

        logger.info(f"Analyzing phenotype {i}/{n_phenotypes}: {pheno_name}")

        # 构建GEMMA LMM命令
        qc_args = ' '.join(config.get_gemma_qc_args())
        notsnp_arg = "-notsnp" if config.gemma.notsnp else ""

        cmd = f"{config.gemma_path} -bfile genotype " \
              f"-k output/kinship.cXX.txt " \
              f"-lmm {config.gemma.lmm_method} " \
              f"-p {pheno_no_header} " \
              f"-n {pheno_col} " \
              f"-c {covariate_file} " \
              f"{qc_args} " \
              f"{notsnp_arg} " \
              f"-o {pheno_name}_lmm"

        log_file = f"gemma_{pheno_name}.log"
        success, error = utils.run_command(cmd, log_file, logger)

        if success:
            assoc_file = os.path.join(gemma_output_dir,
                                      f"{pheno_name}_lmm.assoc.txt")
            if os.path.exists(assoc_file):
                sig_counts = utils.count_significant_snps(
                    assoc_file, [1e-5, 1e-6, 1e-7])
                logger.info(f"  Completed! Significant SNPs: "
                           f"p<1e-5: {sig_counts[0]}, "
                           f"p<1e-6: {sig_counts[1]}, "
                           f"p<1e-7: {sig_counts[2]}")
            else:
                logger.warning(f"  Analysis failed, check log: {log_file}")
        else:
            logger.warning(f"  Analysis failed: {error}")

        logger.info("")

    # 恢复原工作目录
    os.chdir(original_dir)

    return True


def generate_summary_report(config: AnalysisConfig,
                           pheno_names: List[str],
                           logger: logging.Logger) -> bool:
    """
    生成分析汇总报告

    Args:
        config: 分析配置
        pheno_names: 表型名称列表
        logger: 日志对象

    Returns:
        bool: 是否成功
    """
    logger.info("=" * 60)
    logger.info("Generate Analysis Summary")
    logger.info("=" * 60)

    gemma_output_dir = config.get_gemma_output_dir()
    summary_file = os.path.join(config.outdir, "gwas_summary.txt")

    try:
        with open(summary_file, 'w') as f:
            # 写入表头
            f.write("Phenotype_Name,Total_SNPs,")
            f.write("Significant_SNPs_p<1e-5,")
            f.write("Significant_SNPs_p<1e-6,")
            f.write("Significant_SNPs_p<1e-7\n")

            # 统计每个表型的结果
            for i in range(1, len(pheno_names)):  # 跳过第一列（样本ID）
                pheno_name = pheno_names[i]
                assoc_file = os.path.join(gemma_output_dir,
                                          f"{pheno_name}_lmm.assoc.txt")

                if os.path.exists(assoc_file):
                    total_snps = utils.count_lines(assoc_file) - 1  # 减去表头
                    sig_counts = utils.count_significant_snps(
                        assoc_file, [1e-5, 1e-6, 1e-7])

                    f.write(f"{pheno_name},{total_snps},")
                    f.write(f"{sig_counts[0]},{sig_counts[1]},{sig_counts[2]}\n")

        # 显示汇总
        logger.info(f"Summary report: {summary_file}")
        logger.info("")
        with open(summary_file, 'r') as f:
            for line in f:
                logger.info(line.strip())

        return True

    except Exception as e:
        logger.error(f"Failed to generate summary: {str(e)}")
        return False


# ==================== 主函数 ====================

def main():
    """主函数"""
    start_time = time.time()

    # 解析参数
    args = parse_arguments()

    # 设置日志级别
    if args.quiet:
        log_level = logging.ERROR
    elif args.verbose >= 2:
        log_level = logging.DEBUG
    elif args.verbose == 1:
        log_level = logging.INFO
    else:
        log_level = logging.WARNING

    # 初始化日志
    logger = setup_logger(__name__, args.log_file, log_level)

    # 输出程序信息
    logger.info("=" * 60)
    logger.info("Program: GEMMA GWAS Batch Analysis")
    logger.info("Version: 1.0.0")
    logger.info("=" * 60)

    # 创建配置
    config = AnalysisConfig(
        vcf=args.input,
        pheno=args.pheno,
        outdir=args.output_dir,
        n_pca=args.n_pca,
        threads=args.threads,
        gemma_path=args.gemma
    )

    # 设置质控参数
    if args.no_qc:
        config.plink_qc.enable = False
    else:
        if args.maf is not None:
            config.plink_qc.maf = args.maf
        if args.geno is not None:
            config.plink_qc.geno = args.geno
        if args.mind is not None:
            config.plink_qc.mind = args.mind
        if args.hwe is not None:
            config.plink_qc.hwe = args.hwe

    # 设置GEMMA参数
    if args.lmm is not None:
        config.gemma.lmm_method = args.lmm
    if args.gk is not None:
        config.gemma.gk_method = args.gk
    if args.miss_gemma is not None:
        config.gemma.miss = args.miss_gemma
    if args.maf_gemma is not None:
        config.gemma.maf = args.maf_gemma
    config.gemma.notsnp = args.notsnp

    # 验证输入
    logger.info(f"Input VCF: {config.vcf}")
    logger.info(f"Input phenotype: {config.pheno}")
    logger.info(f"Output directory: {config.outdir}")
    logger.info(f"PCA components: {config.n_pca}")
    logger.info(f"LMM method: {config.gemma.lmm_method}")
    logger.info(f"GK method: {config.gemma.gk_method}")

    if config.plink_qc.enable:
        logger.info("Quality control: enabled")
        logger.info(f"  MAF >= {config.plink_qc.maf}")
        logger.info(f"  SNP missing rate <= {config.plink_qc.geno}")
        logger.info(f"  Sample missing rate <= {config.plink_qc.mind}")
        logger.info(f"  HWE p-value >= {config.plink_qc.hwe}")
    else:
        logger.info("Quality control: disabled")

    logger.info(f"GEMMA MAF: {config.gemma.maf}")
    logger.info(f"GEMMA missing rate: {config.gemma.miss}")
    logger.info(f"Threads: {config.threads}")
    logger.info("")

    try:
        # 验证输入文件
        valid, error_msg = config.validate_inputs()
        if not valid:
            logger.critical(f"Input validation failed: {error_msg}")
            sys.exit(1)

        # 检查依赖
        deps_ok, missing = utils.check_dependencies(logger)
        if not deps_ok:
            logger.critical(f"Missing dependencies: {missing}")
            sys.exit(1)

        # 创建输出目录
        os.makedirs(config.outdir, exist_ok=True)

        # 步骤1: 准备表型文件
        success, pheno_files = step1_prepare_phenotype(config, logger)
        if not success:
            sys.exit(1)

        pheno_with_header, pheno_no_header, pheno_names = pheno_files

        # 步骤2: 筛选样本交集（使用bcftools）
        success, pheno_filtered, n_common = step2_filter_common_samples(
            config, pheno_with_header, pheno_no_header, logger
        )
        if not success:
            sys.exit(1)

        # 步骤3: 质量控制（可选）
        if not step3_quality_control(config, logger):
            sys.exit(1)

        # 步骤3b: 修复FAM文件的表型列
        if not step3b_fix_fam_file(config, pheno_filtered, logger):
            sys.exit(1)

        # 步骤4: PCA分析
        success, covariate_file = step4_pca_analysis(
            config, pheno_filtered, logger
        )
        if not success:
            sys.exit(1)

        # 步骤5: 计算亲缘关系矩阵
        if not step5_kinship_matrix(
            config, pheno_filtered, covariate_file, logger
        ):
            sys.exit(1)

        # 步骤6: GWAS分析
        if not step6_gwas_analysis(
            config, pheno_filtered, covariate_file, pheno_names, logger
        ):
            sys.exit(1)

        # 生成汇总报告
        generate_summary_report(config, pheno_names, logger)

        # 输出总结信息
        elapsed_time = time.time() - start_time
        logger.info("=" * 60)
        logger.info("Analysis Summary")
        logger.info("=" * 60)
        logger.info(f"Total runtime: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
        logger.info(f"Output directory: {config.outdir}")
        logger.info(f"Results: {config.get_gemma_output_dir()}/")
        logger.info(f"Summary: {os.path.join(config.outdir, 'gwas_summary.txt')}")
        logger.info("Analysis completed successfully")

    except KeyboardInterrupt:
        logger.warning("Process interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.critical(f"Analysis failed: {str(e)}", exc_info=True)
        sys.exit(1)


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='GEMMA GWAS批量分析 - 使用GEMMA线性混合模型进行GWAS分析|GEMMA GWAS Batch Analysis - Perform GWAS using GEMMA linear mixed model',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i genotype.vcf.gz -p phenotype.txt
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|required arguments')
    required.add_argument('-i', '--input', required=True,
                         help='输入VCF文件(可压缩)|Input VCF file (can be compressed)')
    required.add_argument('-p', '--pheno', required=True,
                         help='表型文件(第一列:样本ID,含表头)|Phenotype file (first column: sample ID, has header)')

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('可选参数|optional arguments')
    optional.add_argument('-o', '--output-dir', default='gemma_results',
                         help='输出目录|Output directory')
    optional.add_argument('--n-pca', type=int, default=10,
                         help='PCA主成分数|Number of PCA components')
    optional.add_argument('--threads', type=int, default=12,
                         help='线程数|Number of threads')
    optional.add_argument('--gemma',
                         default='~/.local/bin/gemma',
                         help='GEMMA程序路径|GEMMA program path')

    # PLINK质控参数|PLINK quality control parameters
    qc_group = parser.add_argument_group('质控参数(PLINK)|quality control parameters (PLINK)')
    qc_group.add_argument('--maf', type=float,
                          help='最小等位基因频率阈值|Minor allele frequency threshold')
    qc_group.add_argument('--geno', type=float,
                          help='SNP缺失率阈值|SNP missing rate threshold')
    qc_group.add_argument('--mind', type=float,
                          help='样本缺失率阈值|Sample missing rate threshold')
    qc_group.add_argument('--hwe', type=float,
                          help='Hardy-Weinberg p值阈值|Hardy-Weinberg p-value threshold')
    qc_group.add_argument('--no-qc', action='store_true',
                          help='跳过PLINK质控|Skip PLINK quality control')

    # GEMMA参数|GEMMA parameters
    gemma_group = parser.add_argument_group('GEMMA参数|GEMMA parameters')
    gemma_group.add_argument('--lmm', type=int, choices=[1, 2, 3, 4],
                             help='LMM检验方法|LMM test method: 1=Wald, 2=LRT, 3=Score, 4=all')
    gemma_group.add_argument('--gk', type=int, choices=[1, 2],
                             help='亲缘关系矩阵方法|Kinship matrix method: 1=centered, 2=standardized')
    gemma_group.add_argument('--miss-gemma', type=float,
                             help='GEMMA缺失率阈值|GEMMA missing rate threshold')
    gemma_group.add_argument('--maf-gemma', type=float,
                             help='GEMMA MAF阈值|GEMMA MAF threshold')
    gemma_group.add_argument('--notsnp', action='store_true',
                             help='不输出每个SNP的估计值(更快)|Do not output estimated values for each SNP (faster)')

    # 日志参数|Logging parameters
    log_group = parser.add_argument_group('日志选项|logging options')
    log_group.add_argument('-v', '--verbose', action='count', default=0,
                          help='详细模式(-v: INFO, -vv: DEBUG)|Verbose mode (-v: INFO, -vv: DEBUG)')
    log_group.add_argument('--quiet', action='store_true',
                          help='静默模式(仅ERROR)|Quiet mode (only ERROR)')
    log_group.add_argument('--log-file',
                          help='日志文件路径|Log file path')

    # 版本信息|Version info
    parser.add_argument('-V', '--version', action='version',
                       version='%(prog)s 1.0.0')

    return parser.parse_args()


if __name__ == '__main__':
    main()
