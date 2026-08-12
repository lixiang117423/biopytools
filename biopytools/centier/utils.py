"""
CentIER工具函数模块|CentIER Utility Functions Module
"""

import logging
import sys
import os
from pathlib import Path
from typing import List
import subprocess


class CentIERLogger:
    """CentIER日志管理器|CentIER Logger Manager"""

    def __init__(self, log_file: Path, log_level: str = "INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.log_file = log_file
        self.log_level = log_level
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 删除已存在的日志文件|Remove existing log file
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, self.log_level.upper(), logging.INFO)

        # stdout handler - INFO级别|stdout handler - INFO level
        # → 超算系统捕获到 .out 文件|→ Captured by job scheduler to .out file
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_formatter = logging.Formatter(log_format, datefmt=date_format)
        stdout_handler.setFormatter(stdout_formatter)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        # → 超算系统捕获到 .err 文件|→ Captured by job scheduler to .err file
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_formatter = logging.Formatter(log_format, datefmt=date_format)
        stderr_handler.setFormatter(stderr_formatter)

        # 文件handler - 所有级别|File handler - all levels
        # → 本地完整日志|→ Local complete log
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter(log_format, datefmt=date_format)
        file_handler.setFormatter(file_formatter)

        # 配置日志|Configure logging
        self.logger = logging.getLogger("CentIER")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False  # 避免重复|Avoid duplicates

        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(stderr_handler)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def check_dependencies(config, logger) -> bool:
    """
    检查依赖工具|Check dependencies

    Args:
        config: CentIERConfig配置对象|CentIERConfig object
        logger: 日志器|Logger

    Returns:
        bool: 是否所有依赖都可用|Whether all dependencies are available
    """
    logger.info("检查依赖工具|Checking dependencies")

    all_ok = True

    # 检查centIER.py脚本|Check centIER.py script
    centier_script = config.get_centier_script_path()
    if not os.path.exists(centier_script):
        logger.error(f"未找到centIER.py脚本|centIER.py script not found: {centier_script}")
        all_ok = False
    else:
        logger.info(f"找到centIER.py|Found centIER.py: {centier_script}")

    # 检查bin目录下的工具|Check tools in bin directory
    bin_path = config.get_bin_path()

    required_tools = [
        ('hmmsearch', 'bin/hmmsearch'),
        ('ltr_finder', 'bin/ltr_finder/ltr_finder'),
        ('trf', 'bin/trf409.linux64'),
        ('REXdb.hmm', 'bin/REXdb.hmm'),
        ('Ty3_gypsy.hmm', 'bin/Ty3_gypsy.hmm')
    ]

    for tool_name, tool_rel_path in required_tools:
        tool_path = os.path.join(config.centier_path, tool_rel_path)
        if os.path.exists(tool_path):
            # 检查可执行权限|Check execute permission
            if tool_rel_path.endswith(('.linux64', 'hmmsearch')) or 'ltr_finder' in tool_rel_path:
                if os.access(tool_path, os.X_OK):
                    logger.info(f"找到且可执行|Found and executable {tool_name}: {tool_path}")
                else:
                    logger.warning(f"找到但不可执行|Found but not executable {tool_name}: {tool_path}")
                    logger.warning(f"请运行|Please run: chmod +x {tool_path}")
                    all_ok = False
            else:
                logger.info(f"找到|Found {tool_name}: {tool_path}")
        else:
            logger.error(f"未找到|{tool_name} not found: {tool_path}")
            all_ok = False

    # 检查外部工具(使用PATH搜索) - gt和LTR_retriever为CentIER必须|
    # Check external tools (use PATH search) - gt and LTR_retriever are required by CentIER
    import shutil as shutil_mod
    external_tools = ['gt', 'LTR_retriever']
    for tool in external_tools:
        tool_path = shutil_mod.which(tool)
        if tool_path:
            logger.info(f"找到|Found {tool}: {tool_path}")
        else:
            # 搜索conda环境中是否有此工具,给出安装建议|Search conda envs for this tool, suggest fix
            conda_envs_dir = os.path.expanduser('~/miniforge3/envs')
            found_envs = []
            if os.path.isdir(conda_envs_dir):
                for env_name in os.listdir(conda_envs_dir):
                    candidate = os.path.join(conda_envs_dir, env_name, 'bin', tool)
                    if os.path.isfile(candidate):
                        found_envs.append(env_name)
            if found_envs:
                logger.error(f"未找到|{tool} not found in PATH (CentIER必须|CentIER requires it)")
                logger.error(f"  但已在以下conda环境中找到|But found in conda envs: {', '.join(found_envs)}")
                logger.error(f"  解决方案|Fix: conda install -n centier -c bioconda genometools")
                logger.error(f"  或|Or: export PATH=$(conda run -n {found_envs[0]} which {tool} | xargs dirname):$PATH")
            else:
                logger.error(f"未找到|{tool} not found in PATH (CentIER必须|CentIER requires it)")
                logger.error(f"  解决方案|Fix: conda install -n centier -c bioconda genometools")
            all_ok = False

    # Hi-C FASTQ 自动模式额外校验|Extra checks for Hi-C FASTQ auto mode
    if getattr(config, 'fastq_r1', None):
        from ..common.paths import get_tool_path as _get_tool_path
        hicpro_path = _get_tool_path(
            'hicpro',
            '~/software/HiC-Pro_v3.1.0/HiC-Pro_3.1.0/bin/HiC-Pro',
            'HICPRO_PATH'
        )
        if os.path.exists(hicpro_path):
            logger.info(f"找到 HiC-Pro|Found HiC-Pro: {hicpro_path}")
        else:
            logger.error(f"未找到 HiC-Pro|HiC-Pro not found: {hicpro_path}")
            logger.error("  Hi-C FASTQ 自动模式需要 HiC-Pro|Hi-C FASTQ auto mode requires HiC-Pro")
            all_ok = False

        # bowtie2-build 在 HiC-Pro env 内|bowtie2-build lives in HiC-Pro env
        bowtie2_build_path = _get_tool_path(
            'bowtie2_build',
            '~/miniforge3/envs/HiC-Pro_v3.1.0/bin/bowtie2-build',
            'BOWTIE2_BUILD_PATH'
        )
        if os.path.exists(bowtie2_build_path):
            logger.info(f"找到 bowtie2-build|Found bowtie2-build: {bowtie2_build_path}")
        else:
            # 缺失只警告,不硬失败(HiC-Pro 可能用预建索引)|Warn only, don't hard-fail
            logger.warning(f"bowtie2-build 未找到(自动建索引可能失败)|"
                           f"bowtie2-build not found (auto index build may fail): {bowtie2_build_path}")

    # 检查Python包依赖|Check Python package dependencies
    import importlib.util
    python_packages = ['pyfastx', 'numpy', 'pandas', 'scipy']
    for pkg in python_packages:
        spec = importlib.util.find_spec(pkg)
        if spec is not None:
            logger.info(f"找到Python包|Found Python package {pkg}: {spec.origin}")
        else:
            logger.warning(f"未找到Python包|Python package {pkg} not found (CentIER需要|CentIER requires it)")

    if all_ok:
        logger.info("依赖检查完成|Dependency check completed")
    else:
        logger.error("依赖检查失败|Dependency check failed")

    return all_ok


def run_command(cmd: List[str], logger, check: bool = True) -> subprocess.CompletedProcess:
    """
    运行命令|Run command

    Args:
        cmd: 命令列表|Command list
        logger: 日志器|Logger
        check: 是否检查返回码|Whether to check return code

    Returns:
        subprocess.CompletedProcess: 命令执行结果|Command execution result
    """
    logger.info(f"命令|Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=check
        )
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败|Command execution failed: {e}")
        logger.error(f"标准错误|Stderr: {e.stderr}")
        raise


def format_number(num: int) -> str:
    """格式化数字|Format number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)
