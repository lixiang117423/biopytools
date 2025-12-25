"""
🔧 GWAS分析工具函数模块 | GWAS Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional


class GWASLogger:
    """GWAS分析日志管理器 | GWAS Analysis Logger Manager"""

    def __init__(self, output_prefix: str, log_name: str = "gwas_analysis.log"):
        self.output_prefix = output_prefix
        self.log_file = f"{output_prefix}.log"
        self.setup_logging()

    def setup_logging(self):
        """设置日志 | Setup logging"""
        # 如果日志文件存在则删除 | Remove log file if exists
        if Path(self.log_file).exists():
            Path(self.log_file).unlink()

        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout),
            ],
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger


class CommandRunner:
    """命令执行器 | Command Runner"""

    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir.resolve()

    def run(self, cmd: str, description: str = "", timeout: int = 3600) -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"🚀 执行步骤 | Executing step: {description}")

        self.logger.info(f"📋 命令 | Command: {cmd}")
        self.logger.info(f"📁 工作目录 | Working directory: {self.working_dir}")

        start_time = time.time()

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.working_dir,
                timeout=timeout,
            )

            end_time = time.time()
            duration = end_time - start_time
            self.logger.info(
                f"✅ 命令执行成功 | Command executed successfully: {description} (耗时 | Duration: {duration:.2f}s)"
            )

            if result.stdout:
                self.logger.debug(f"📤 标准输出 | Stdout: {result.stdout}")

            return True

        except subprocess.CalledProcessError as e:
            end_time = time.time()
            duration = end_time - start_time
            self.logger.error(
                f"❌ 命令执行失败 | Command execution failed: {description} (耗时 | Duration: {duration:.2f}s)"
            )
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"💥 错误信息 | Error message: {e.stderr}")
            if e.stdout:
                self.logger.error(f"📤 标准输出 | Stdout: {e.stdout}")
            return False

        except subprocess.TimeoutExpired:
            self.logger.error(
                f"⏰ 命令执行超时 | Command execution timeout: {description} (超时 | Timeout: {timeout}s)"
            )
            return False


def check_dependencies(config, logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("🔍 检查依赖软件 | Checking dependencies")

    dependencies = [
        (config.plink_path, "PLINK"),
        (config.bcftools_path, "BCFtools"),
        (config.admixture_path, "ADMIXTURE"),
        (config.emmax_path, "EMMAX"),
    ]

    missing_deps = []

    for cmd, name in dependencies:
        try:
            result = subprocess.run(
                [cmd, "--help"], capture_output=True, text=True, timeout=10
            )
            if (
                result.returncode == 0
                or "usage" in result.stdout.lower()
                or "usage" in result.stderr.lower()
            ):
                logger.info(f"✅ {name} 可用 | available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)

    # 检查Python依赖 | Check Python dependencies
    python_deps = [
        ("pandas", "pandas"),
        ("numpy", "numpy"),
        ("matplotlib", "matplotlib"),
        ("seaborn", "seaborn"),
        ("scipy", "scipy"),
    ]

    python_missing = []
    for dep_name, import_name in python_deps:
        try:
            __import__(import_name)
            logger.info(
                f"✅ Python包 {dep_name} 可用 | Python package {dep_name} available"
            )
        except ImportError:
            python_missing.append(dep_name)

    if missing_deps:
        error_msg = f"❌ 缺少外部依赖软件 | Missing external dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)

    if python_missing:
        error_msg = f"❌ 缺少Python依赖包 | Missing Python dependencies: {', '.join(python_missing)}"
        logger.error(error_msg)

        # 提供安装建议 | Provide installation suggestions
        logger.error("📦 安装建议 | Installation suggestions:")
        for dep in python_missing:
            logger.error(f"  pip install {dep}")

        raise RuntimeError(error_msg)

    logger.info("✅ 所有依赖检查通过 | All dependencies check passed")
    return True


def get_vcf_stats(vcf_file: str, logger) -> Dict[str, Any]:
    """获取VCF文件统计信息 | Get VCF file statistics"""
    try:
        stats = {}

        # 统计样本数量 | Count samples
        with open(vcf_file, "r") as f:
            for line in f:
                if line.startswith("#CHROM"):
                    headers = line.strip().split("\t")
                    # 前9列是固定的，从第10列开始是样本 | First 9 columns are fixed, samples start from column 10
                    stats["samples"] = len(headers) - 9
                    break

        # 统计变异数量 | Count variants
        variant_count = 0
        with open(vcf_file, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    variant_count += 1
        stats["variants"] = variant_count

        # 统计染色体信息 | Count chromosome information
        chromosomes = set()
        with open(vcf_file, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    parts = line.strip().split("\t")
                    if parts:
                        chromosomes.add(parts[0])
        stats["chromosomes"] = sorted(list(chromosomes))

        return stats

    except Exception as e:
        logger.warning(f"⚠️ 无法获取VCF统计信息 | Cannot get VCF statistics: {e}")
        return {}


def get_phenotype_stats(phenotype_file: str, logger) -> Dict[str, Any]:
    """获取表型文件统计信息 | Get phenotype file statistics"""
    try:
        stats = {}

        # 读取表型文件 | Read phenotype file
        df = pd.read_csv(
            phenotype_file, sep="\t", header=None, names=["sample_id", "phenotype"]
        )

        stats["samples"] = len(df)
        stats["phenotype_mean"] = df["phenotype"].mean()
        stats["phenotype_std"] = df["phenotype"].std()
        stats["phenotype_min"] = df["phenotype"].min()
        stats["phenotype_max"] = df["phenotype"].max()
        stats["missing_values"] = df["phenotype"].isna().sum()

        return stats

    except Exception as e:
        logger.warning(f"⚠️ 无法获取表型统计信息 | Cannot get phenotype statistics: {e}")
        return {}


def create_temp_dir(base_dir: str, prefix: str = "gwas_temp") -> Path:
    """创建临时目录 | Create temporary directory"""
    temp_dir = Path(base_dir) / f"{prefix}_{int(time.time())}"
    temp_dir.mkdir(parents=True, exist_ok=True)
    return temp_dir


def cleanup_temp_files(temp_dirs: List[Path], logger):
    """清理临时文件 | Cleanup temporary files"""
    for temp_dir in temp_dirs:
        try:
            if temp_dir.exists():
                import shutil

                shutil.rmtree(temp_dir)
                logger.info(
                    f"🗑️ 已清理临时目录 | Cleaned temporary directory: {temp_dir}"
                )
        except Exception as e:
            logger.warning(
                f"⚠️ 无法清理临时目录 | Cannot clean temporary directory: {temp_dir}, 错误: {e}"
            )
