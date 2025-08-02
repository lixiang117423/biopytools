"""
工具模块 | Utilities Module
包含日志、命令执行、依赖检查等通用功能
"""

import glob
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple


class MethylationLogger:
    """甲基化分析日志管理器 | Methylation Analysis Logger Manager"""

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.log_dir = self.output_dir / "logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)

        # 设置日志文件路径 | Set log file path
        from datetime import datetime

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.log_file = self.log_dir / f"methylation_pipeline_{timestamp}.log"

        # 配置日志 | Configure logging
        self.logger = logging.getLogger("methylation_pipeline")
        self.logger.setLevel(logging.INFO)

        # 避免重复添加handler | Avoid duplicate handlers
        if not self.logger.handlers:
            # 文件处理器 | File handler
            file_handler = logging.FileHandler(self.log_file, encoding="utf-8")
            file_handler.setLevel(logging.INFO)

            # 控制台处理器 | Console handler
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setLevel(logging.INFO)

            # 设置格式 | Set format
            formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s")
            file_handler.setFormatter(formatter)
            console_handler.setFormatter(formatter)

            self.logger.addHandler(file_handler)
            self.logger.addHandler(console_handler)

    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器 | Command Runner"""

    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = Path(working_dir)
        self.working_dir.mkdir(parents=True, exist_ok=True)

    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"🔧 执行步骤 | Executing step: {description}")

        self.logger.info(f"命令 | Command: {cmd}")
        self.logger.info(f"工作目录 | Working directory: {self.working_dir}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir,
            )

            self.logger.info(
                f"✅ 命令执行成功 | Command executed successfully: {description}"
            )

            if result.stdout:
                self.logger.debug(f"标准输出 | Stdout: {result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(
                f"❌ 命令执行失败 | Command execution failed: {description}"
            )
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            self.logger.error(f"错误信息 | Error message: {e.stderr}")
            self.logger.error(f"标准输出 | Stdout: {e.stdout}")
            return False


def check_dependencies(config, logger):
    """🔍 检查依赖软件 | Check dependencies"""
    logger.info("🔍 检查依赖软件 | Checking dependencies")

    # 必需工具 | Required tools
    required_tools = [
        (config.fastp_path, "fastp"),
        (config.bismark_path, "bismark"),
        (config.bismark_genome_preparation_path, "bismark_genome_preparation"),
        (config.bowtie2_path, "bowtie2"),
        (config.deduplicate_bismark_path, "deduplicate_bismark"),
        (config.bismark_methylation_extractor_path, "bismark_methylation_extractor"),
        (config.bismark2report_path, "bismark2report"),
        (config.bismark2summary_path, "bismark2summary"),
    ]

    # 可选工具 | Optional tools
    optional_tools = [
        (config.makeblastdb_path, "makeblastdb"),
        (config.blastn_path, "blastn"),
        (config.multiqc_path, "multiqc"),
        (config.bedtools_path, "bedtools"),
        (config.samtools_path, "samtools"),
    ]

    missing_required = []
    missing_optional = []
    python_deps_missing = False

    # 检查必需工具 | Check required tools
    for cmd, name in required_tools:
        if not shutil.which(cmd):
            missing_required.append(name)
            logger.error(f"❌ 错误: {name} 未安装或不在PATH中")
            logger.error(f"   配置的路径: {cmd}")
            logger.error(f"   安装建议: conda install -c bioconda {name.lower()}")
        else:
            logger.info(f"✅ {name}: {shutil.which(cmd)}")

    # 检查可选工具 | Check optional tools
    for cmd, name in optional_tools:
        if not shutil.which(cmd):
            missing_optional.append(name)
            logger.warning(f"⚠️  {name} 未安装，部分功能可能受限")
            logger.warning(f"    配置的路径: {cmd}")
            logger.warning(f"    安装建议: conda install -c bioconda {name.lower()}")
        else:
            logger.info(f"✅ {name}: {shutil.which(cmd)}")

    # 检查R环境 | Check R environment
    if config.enhanced_mode:
        if not os.path.exists(config.r_executable):
            missing_required.append("R")
            logger.error(f"❌ 错误: R 未找到，路径: {config.r_executable}")
            logger.error("   请检查R安装或修改--r-executable参数")
        else:
            logger.info(f"✅ R: {config.r_executable}")

            # 检查R包 | Check R packages
            r_packages = [
                "methylKit",
                "GenomicRanges",
                "rtracklayer",
                "ggplot2",
                "dplyr",
            ]
            r_lib_path = (
                "/share/org/YZWL/yzwl_lixg/miniforge3/envs/methylkit/lib/R/library"
            )

            for pkg in r_packages:
                try:
                    result = subprocess.run(
                        [
                            config.r_executable,
                            "--slave",
                            "-e",
                            f".libPaths('{r_lib_path}'); library({pkg}); cat('{pkg}_OK')",
                        ],
                        capture_output=True,
                        text=True,
                        timeout=30,
                    )

                    if f"{pkg}_OK" in result.stdout:
                        logger.info(f"✅ {pkg} R package 已安装")
                    else:
                        logger.warning(f"⚠️  {pkg} R package 检测异常，部分功能可能受限")
                        logger.warning(
                            f"    安装建议: install.packages('{pkg}') 或 BiocManager::install('{pkg}')"
                        )

                except (subprocess.TimeoutExpired, FileNotFoundError):
                    logger.warning(f"⚠️  {pkg} R package 检测失败")

    # 检查Python依赖包 | Check Python dependencies
    python_deps_missing = False
    try:
        import json

        logger.info("✅ Python json 模块可用")
    except ImportError:
        logger.warning("⚠️  Python json 模块不可用，JSON样品分组文件功能受限")
        python_deps_missing = True

    try:
        import numpy as np
        import pandas as pd

        logger.info("✅ Python pandas, numpy 模块可用")
    except ImportError:
        logger.warning("⚠️  Python pandas/numpy 模块不可用，部分增强分析功能受限")
        logger.warning("    安装建议: pip install pandas numpy")
        if config.enhanced_mode:
            python_deps_missing = True

    # 总结检查结果 | Summarize check results
    if missing_required:
        error_msg = f"缺少必需软件 | Missing required dependencies: {', '.join(missing_required)}"
        logger.error(error_msg)
        logger.error(
            "请安装缺少的软件后重新运行 | Please install missing software and try again"
        )
        logger.error("常用安装命令 | Common installation commands:")
        logger.error("  conda install -c bioconda fastp bismark bowtie2")
        logger.error("  conda install -c conda-forge r-base")
        raise RuntimeError(error_msg)

    if python_deps_missing and config.enhanced_mode:
        error_msg = (
            "增强分析模式需要Python依赖包 | Enhanced mode requires Python dependencies"
        )
        logger.error(error_msg)
        logger.error("安装命令: pip install pandas numpy")
        raise RuntimeError(error_msg)

    if missing_optional:
        logger.warning(
            f"缺少可选软件 | Missing optional dependencies: {', '.join(missing_optional)}"
        )
        logger.warning("这些软件缺失可能影响部分高级功能")
        logger.warning("建议安装命令 | Recommended installation:")
        logger.warning("  conda install -c bioconda blast bedtools samtools multiqc")

    logger.info("✅ 依赖软件检查完成")
    return True


def clean_sample_name(sample_name: str) -> str:
    """清理样品名称（去掉前缀）| Clean sample name (remove prefix)"""
    sample_name = sample_name.replace("FZYM412_", "")
    return sample_name


def check_completed_steps(config, logger) -> dict:
    """🔍 检测已完成的步骤 | Check completed steps"""
    logger.info("🔍 检测已完成的步骤...")

    status = {
        "need_rename": False,
        "fastp_completed": True,
        "index_exists": False,
        "completed_extraction": [],
        "expected_samples": [],
    }

    # 检查是否需要文件重命名
    if os.path.exists(config.raw_dir):
        for file in os.listdir(config.raw_dir):
            if "?" in file or file.startswith("FZYM412_"):
                status["need_rename"] = True
                break

    # 获取预期样品列表 - 改进的样品名提取逻辑
    if os.path.exists(config.raw_dir):
        sample_files = {}

        for filename in os.listdir(config.raw_dir):
            if filename.endswith(".fq.gz"):
                if "_1.fq.gz" in filename:
                    base_name = filename.replace("_1.fq.gz", "")
                    if base_name not in sample_files:
                        sample_files[base_name] = {}
                    sample_files[base_name]["R1"] = filename
                elif "_2.fq.gz" in filename:
                    base_name = filename.replace("_2.fq.gz", "")
                    if base_name not in sample_files:
                        sample_files[base_name] = {}
                    sample_files[base_name]["R2"] = filename

        # 只保留有完整配对文件的样品，并清理样品名
        raw_samples = []
        for base_name, files in sample_files.items():
            if "R1" in files and "R2" in files:
                clean_sample_name = base_name.replace("FZYM412_", "")
                raw_samples.append(clean_sample_name)

        status["expected_samples"] = raw_samples

    # 检查fastp是否已完成
    for sample in status["expected_samples"]:
        clean_r1 = os.path.join(config.clean_dir, f"{sample}_R1_clean.fq.gz")
        clean_r2 = os.path.join(config.clean_dir, f"{sample}_R2_clean.fq.gz")

        if not (os.path.exists(clean_r1) and os.path.exists(clean_r2)):
            status["fastp_completed"] = False
            break

    # 检查Bismark索引
    bismark_index_dir = os.path.join(config.mapping_dir, "bismark_index")
    index_file = os.path.join(
        bismark_index_dir,
        "Bisulfite_Genome",
        "CT_conversion",
        "genome_mfa.CT_conversion.fa",
    )
    status["index_exists"] = os.path.exists(index_file)

    # 检查甲基化提取
    for sample in status["expected_samples"]:
        sample_dir = os.path.join(config.mapping_dir, "bismark_results", sample)
        if os.path.exists(sample_dir):
            for file in os.listdir(sample_dir):
                if file.endswith("CX_report.txt"):
                    status["completed_extraction"].append(sample)
                    break

    # 输出检测结果
    logger.info("🎯 步骤完成情况检测结果:")
    logger.info("==================================")
    logger.info(
        f"📝 步骤1 (文件重命名): {'❌ 需要执行' if status['need_rename'] else '✅ 已完成或无需执行'}"
    )
    logger.info(
        f"🧹 步骤2 (fastp质控): {'✅ 已完成' if status['fastp_completed'] else '❌ 需要执行'}"
    )
    logger.info(
        f"🏗️ 步骤3 (Bismark索引): {'✅ 已完成' if status['index_exists'] else '❌ 需要执行'}"
    )
    logger.info(
        f"🧬 步骤4 (甲基化分析): {len(status['completed_extraction'])}/{len(status['expected_samples'])} 个样品已完成"
    )

    if status["expected_samples"]:
        logger.info(f"📋 检测到的样品: {', '.join(status['expected_samples'])}")

    return status
