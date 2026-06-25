"""
VG工具函数|VG Utility Functions
"""

import os
import subprocess
import logging
from pathlib import Path
from datetime import datetime
from typing import Tuple, Optional


class VGLogger:
    """VG日志类|VG Logger Class"""

    def __init__(self, log_file: Path, log_level: str = "INFO"):
        """
        初始化日志|Initialize logger

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR)
        """
        import sys

        self.log_file = Path(log_file)
        self.log_level = getattr(logging, log_level.upper(), logging.INFO)

        # 确保日志目录存在|Ensure log directory exists
        self.log_file.parent.mkdir(parents=True, exist_ok=True)

        # 创建logger|Create logger
        self.logger = logging.getLogger("VG")
        self.logger.setLevel(logging.DEBUG)

        # 清除已有的处理器|Clear existing handlers
        self.logger.handlers = []
        self.logger.propagate = False

        # 格式化器|Formatter
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file, mode='w', encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    def get_logger(self) -> logging.Logger:
        """获取logger对象|Get logger object"""
        return self.logger


def validate_vg_environment(vg_env: str, logger: Optional[logging.Logger] = None) -> bool:
    """
    验证VG环境|Validate VG environment

    Args:
        vg_env: conda环境名称或路径|conda env name or path
        logger: 日志对象|Logger object

    Returns:
        是否可用|Whether available
    """
    if logger:
        logger.info(" 检查VG环境|Checking VG environment")

    # 构建conda run命令|Build conda run command
    cmd = f"conda run -n {vg_env} --no-capture-output vg --version"

    if logger:
        logger.info(f"   命令|Command: {cmd}")

    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode == 0:
            version_info = result.stdout.strip()
            if logger:
                logger.info(f"   VG可用|VG available: {version_info}")
            return True
        else:
            if logger:
                logger.error(f"   VG不可用|VG not available: {result.stderr}")
            return False

    except Exception as e:
        if logger:
            logger.error(f"   VG环境检查失败|VG environment check failed: {e}")
        return False


def run_vg_command(
    vg_env: str,
    command: str,
    logger: Optional[logging.Logger] = None,
    cwd: Optional[str] = None
) -> Tuple[bool, str, str]:
    """
    运行VG命令|Run VG command

    Args:
        vg_env: conda环境名称|conda env name
        command: VG命令（不含vg前缀）|VG command (without 'vg' prefix)
        logger: 日志对象|Logger object
        cwd: 工作目录|Working directory

    Returns:
        (是否成功, stdout, stderr)|(success, stdout, stderr)
    """
    # 构建完整命令|Build full command
    full_command = f"conda run -n {vg_env} --no-capture-output {command}"

    if logger:
        logger.info(f"   执行命令|Executing command: {full_command}")

    try:
        result = subprocess.run(
            full_command,
            shell=True,
            capture_output=True,
            text=True,
            cwd=cwd,
            timeout=None  # VG可能需要很长时间|VG may take a long time
        )

        success = result.returncode == 0

        if not success and logger:
            logger.error(f"   命令失败，返回码|Command failed, return code: {result.returncode}")
            if result.stderr:
                logger.error(f"   错误信息|Error message: {result.stderr[:500]}")

        return success, result.stdout, result.stderr

    except subprocess.TimeoutExpired:
        if logger:
            logger.error("   命令执行超时|Command execution timeout")
        return False, "", "Timeout"
    except Exception as e:
        if logger:
            logger.error(f"   命令执行异常|Command execution exception: {e}")
        return False, "", str(e)


def validate_reference_file(reference: str, logger: Optional[logging.Logger] = None) -> bool:
    """
    验证参考基因组文件|Validate reference genome file

    Args:
        reference: 参考基因组文件路径|Reference genome file path
        logger: 日志对象|Logger object

    Returns:
        是否有效|Whether valid
    """
    if not Path(reference).exists():
        if logger:
            logger.error(f"参考基因组文件不存在|Reference file does not exist: {reference}")
        return False

    # 检查文件扩展名|Check file extension
    if not reference.endswith(('.fa', '.fasta', '.fna')):
        if logger:
            logger.warning(f"参考基因组文件扩展名可能不正确|Reference file extension may be incorrect: {reference}")

    if logger:
        logger.info(f"   参考基因组文件验证通过|Reference file validated: {reference}")

    return True


def validate_vcf_file(vcf: str, logger: Optional[logging.Logger] = None) -> bool:
    """
    验证VCF文件|Validate VCF file

    Args:
        vcf: VCF文件路径|VCF file path
        logger: 日志对象|Logger object

    Returns:
        是否有效|Whether valid
    """
    if not Path(vcf).exists():
        if logger:
            logger.error(f"VCF文件不存在|VCF file does not exist: {vcf}")
        return False

    # 检查是否压缩和索引|Check if compressed and indexed
    if vcf.endswith('.gz'):
        index_file = f"{vcf}.tbi"
        if not Path(index_file).exists():
            if logger:
                logger.warning(f"VCF索引文件不存在|VCF index file does not exist: {index_file}")
        else:
            if logger:
                logger.info(f"   VCF文件已索引|VCF file indexed: {vcf}")
    else:
        if logger:
            logger.warning(f"VCF文件未压缩，建议使用bgzip压缩|VCF file not compressed, bgzip recommended: {vcf}")

    if logger:
        logger.info(f"   VCF文件验证通过|VCF file validated: {vcf}")

    return True


def check_output_files(output_files: list) -> bool:
    """
    检查输出文件是否存在|Check if output files exist

    Args:
        output_files: 输出文件列表|List of output files

    Returns:
        是否所有文件都存在|Whether all files exist
    """
    for output_file in output_files:
        if not Path(output_file).exists():
            return False
    return True


def get_file_size_mb(file_path: str) -> float:
    """
    获取文件大小（MB）|Get file size in MB

    Args:
        file_path: 文件路径|File path

    Returns:
        文件大小（MB）|File size in MB
    """
    if not Path(file_path).exists():
        return 0.0

    size_bytes = Path(file_path).stat().st_size
    return size_bytes / (1024 * 1024)
