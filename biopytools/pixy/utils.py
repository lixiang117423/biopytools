"""
Pixy工具模块|Pixy Utilities Module
"""

import logging
import sys
from pathlib import Path
from typing import Optional
import subprocess


class PixyLogger:
    """Pixy日志管理类|Pixy Logger Management Class"""

    def __init__(self, log_file: Path):
        """初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
        """
        self.log_file = log_file
        self.logger = self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        """设置日志器|Setup logger

        Returns:
            配置好的日志器|Configured logger
        """
        logger = logging.getLogger('Pixy')
        logger.setLevel(logging.DEBUG)

        # 清除已有的处理器|Clear existing handlers
        logger.handlers.clear()

        # 文件处理器|File handler
        fh = logging.FileHandler(self.log_file, encoding='utf-8')
        fh.setLevel(logging.DEBUG)

        # 控制台处理器|Console handler
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.INFO)

        # 格式化器|Formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)

        logger.addHandler(fh)
        logger.addHandler(ch)

        return logger

    def get_logger(self) -> logging.Logger:
        """获取日志器|Get logger

        Returns:
            日志器对象|Logger object
        """
        return self.logger


class PixyChecker:
    """Pixy环境检查类|Pixy Environment Checker"""

    def __init__(self, logger, pixy_path: str, conda_env: str):
        """初始化检查器|Initialize checker

        Args:
            logger: 日志器|Logger
            pixy_path: pixy可执行文件路径|pixy executable path
            conda_env: conda环境路径|conda environment path
        """
        self.logger = logger
        self.pixy_path = pixy_path
        self.conda_env = conda_env

    def check_pixy(self) -> bool:
        """检查pixy是否可用|Check if pixy is available

        Returns:
            bool: pixy是否可用|Whether pixy is available
        """
        self.logger.info(f"检查pixy环境|Checking pixy environment")

        # 首先检查conda环境是否存在|First check if conda environment exists
        if not Path(self.conda_env).exists():
            self.logger.error(f"Conda环境不存在|Conda environment does not exist: {self.conda_env}")
            return False

        # 从conda环境路径提取环境名|Extract environment name from conda environment path
        env_name = Path(self.conda_env).name

        # 检查pixy命令|Check pixy command
        try:
            # 尝试在conda环境中运行pixy --version|Try to run pixy --version in conda environment
            result = subprocess.run(
                f"conda run -n {env_name} pixy --version",
                shell=True,
                capture_output=True,
                text=True,
                timeout=30
            )

            if result.returncode == 0:
                version = result.stdout.strip()
                self.logger.info(f"Pixy版本|Pixy version: {version}")
                return True
            else:
                self.logger.error(f"Pixy不可用|Pixy is not available in conda environment: {self.conda_env}")
                return False

        except subprocess.TimeoutExpired:
            self.logger.error("Pixy版本检查超时|Pixy version check timeout")
            return False
        except Exception as e:
            self.logger.error(f"Pixy检查失败|Pixy check failed: {e}")
            return False

    def check_dependencies(self) -> bool:
        """检查pixy依赖是否满足|Check if pixy dependencies are met

        Returns:
            bool: 依赖是否满足|Whether dependencies are met
        """
        self.logger.info("检查pixy依赖|Checking pixy dependencies")

        # pixy依赖bcftools和tabix|pixy depends on bcftools and tabix
        dependencies = ['bcftools', 'tabix']

        # 从conda环境路径提取环境名|Extract environment name from conda environment path
        env_name = Path(self.conda_env).name

        for dep in dependencies:
            try:
                result = subprocess.run(
                    f"conda run -n {env_name} {dep} --version",
                    shell=True,
                    capture_output=True,
                    text=True,
                    timeout=10
                )

                if result.returncode == 0:
                    self.logger.info(f"  {dep} 可用|{dep} is available")
                else:
                    self.logger.error(f"  {dep} 不可用|{dep} is not available")
                    return False

            except subprocess.TimeoutExpired:
                self.logger.error(f"  {dep} 检查超时|{dep} check timeout")
                return False
            except Exception as e:
                self.logger.error(f"  {dep} 检查失败|{dep} check failed: {e}")
                return False

        return True


def parse_population_file(pop_file: Path) -> dict:
    """解析群体文件|Parse population file

    Args:
        pop_file: 群体文件路径|Population file path

    Returns:
        dict: 样本到群体的映射字典|Sample to population mapping dictionary
    """
    pop_map = {}
    with open(pop_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            sample_id = parts[0]
            population = parts[1]
            pop_map[sample_id] = population

    return pop_map


def check_vcf_has_invariant_sites(vcf_file: Path, conda_env: str, logger) -> bool:
    """检查VCF是否包含不变位点|Check if VCF contains invariant sites

    Args:
        vcf_file: VCF文件路径|VCF file path
        conda_env: conda环境路径|conda environment path
        logger: 日志器|Logger

    Returns:
        bool: 是否包含不变位点|Whether contains invariant sites
    """
    logger.info("检查VCF是否包含不变位点|Checking if VCF contains invariant sites")

    # 从conda环境路径提取环境名|Extract environment name from conda environment path
    env_name = Path(conda_env).name

    try:
        # 使用bcftools检查前100个位点是否有ALT="."的位点（减少检查量避免超时）
        # Use bcftools to check if there are sites with ALT="." in first 100 sites
        result = subprocess.run(
            f"conda run -n {env_name} bcftools view -H {vcf_file} | head -n 100 | grep -c '\\t\\.\\t'",
            shell=True,
            capture_output=True,
            text=True,
            timeout=120  # 增加超时时间到120秒|Increase timeout to 120s
        )

        # grep返回0表示找到了不变位点，返回1表示没找到
        # grep returns 0 if found, 1 if not found
        has_invariant = result.returncode == 0

        if has_invariant:
            count = int(result.stdout.strip())
            logger.info(f"  VCF包含不变位点（前100个位点中发现{count}个）|VCF contains invariant sites (found {count} in first 100 sites)")
        else:
            logger.warning("  ⚠ VCF不包含不变位点（ALT=\".\"），将自动绕过检查|VCF does not contain invariant sites (ALT=\".\"), will bypass check automatically")
            logger.warning("  ⚠ 注意：这会导致pi和dxy的统计结果不准确|Note: This may cause inaccurate pi and dxy results")

        return has_invariant

    except subprocess.TimeoutExpired:
        # 超时说明VCF很大且可能不包含不变位点，假设不包含
        # Timeout means VCF is large and likely doesn't contain invariant sites, assume no invariant sites
        logger.warning("  ⚠ VCF不变位点检查超时，假设不包含不变位点|VCF invariant sites check timeout, assuming no invariant sites")
        logger.warning("  ⚠ 将自动绕过检查|Will bypass check automatically")
        return False
    except Exception as e:
        # 出错时保守起见，假设不包含不变位点
        # On error, be conservative and assume no invariant sites
        logger.warning(f"  ⚠ VCF不变位点检查失败: {e}，假设不包含不变位点|VCF invariant sites check failed: {e}, assuming no invariant sites")
        return False


def get_population_summary(pop_map: dict) -> dict:
    """获取群体统计信息|Get population summary

    Args:
        pop_map: 样本到群体的映射字典|Sample to population mapping dictionary

    Returns:
        dict: 群体统计信息|Population summary information
    """
    pop_counts = {}
    for sample, pop in pop_map.items():
        pop_counts[pop] = pop_counts.get(pop, 0) + 1

    return {
        'total_samples': len(pop_map),
        'total_populations': len(pop_counts),
        'population_counts': pop_counts
    }
