"""
PGGB工具函数|PGGB Utility Functions
"""

import os
import subprocess
import logging
import sys
from pathlib import Path
from typing import Optional, List

# conda包装统一走公共层(§13): 同源conda绝对路径 + run -p <环境前缀>, 严禁裸调conda
from ..common.conda_runner import build_conda_command


class PGGBLogger:
    """PGGB日志管理器|PGGB Logger Manager"""

    def __init__(self, log_file, log_level: str = "INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level
        """
        self.log_file = Path(log_file)
        self.log_level = getattr(logging, log_level.upper(), logging.INFO)

        # 确保日志目录存在|Ensure log directory exists
        self.log_file.parent.mkdir(parents=True, exist_ok=True)

        # 创建logger|Create logger
        self.logger = logging.getLogger("PGGB")
        self.logger.setLevel(logging.DEBUG)
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


def check_pggb_dependencies(config, logger: logging.Logger) -> bool:
    """
    检查PGGB依赖工具是否可用|Check PGGB dependency tools availability

    Args:
        config: PGGBConfig配置对象|PGGBConfig configuration object
        logger: 日志对象|Logger object

    Returns:
        是否全部可用|Whether all dependencies are available
    """
    tools = ['pggb']

    # 如果需要VCF输出，检查额外工具
    if config.vcf_spec:
        tools.extend(['vg', 'bcftools'])

    all_ok = True
    for tool in tools:
        cmd = build_conda_command(tool, ['--version'])
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            version_line = result.stdout.strip()
            if result.returncode == 0 and version_line:
                logger.info(f"  {tool}: {version_line.split(chr(10))[0]}")
            else:
                # pggb的--version输出到stderr
                version_line = result.stderr.strip()
                if version_line:
                    logger.info(f"  {tool}: {version_line.split(chr(10))[0]}")
                else:
                    logger.warning(f"  {tool}: 无法获取版本信息|Cannot get version info")
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.error(f"  {tool}: 不可用|Not available - {e}")
            all_ok = False

    return all_ok


def ensure_fasta_index(fasta_path: str, logger: logging.Logger) -> bool:
    """
    确保FASTA索引文件存在|Ensure FASTA index file exists

    Args:
        fasta_path: FASTA文件路径|FASTA file path
        logger: 日志对象|Logger object

    Returns:
        是否成功|Whether successful
    """
    fai_path = fasta_path + '.fai'

    if Path(fai_path).exists():
        logger.info(f"  FASTA索引已存在|FASTA index already exists: {fai_path}")
        return True

    logger.info(f"  生成FASTA索引|Generating FASTA index...")
    cmd = build_conda_command('samtools', ['faidx', fasta_path])
    logger.info(f"  命令|Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=300)
        if Path(fai_path).exists():
            logger.info(f"  索引生成成功|Index generated: {fai_path}")
            return True
        else:
            logger.error(f"  索引生成失败|Index generation failed")
            return False
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        logger.error(f"  索引生成失败|Index generation failed: {e}")
        return False


def format_number(num: int) -> str:
    """
    格式化数字为大单位显示|Format number with large unit display

    Args:
        num: 数字|Number

    Returns:
        格式化后的字符串|Formatted string
    """
    if num >= 1_000_000_000:
        return f"{num / 1_000_000_000:.2f}G"
    elif num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.1f}K"
    return str(num)


def format_size(size_bytes: int) -> str:
    """
    将字节数格式化为人类可读的大小|Format bytes to human-readable size

    Args:
        size_bytes: 字节数|Number of bytes

    Returns:
        格式化后的字符串|Formatted string
    """
    if size_bytes >= 1024 ** 3:
        size_gb = size_bytes / (1024 ** 3)
        return f"{size_gb:.0f}G" if size_gb >= 10 else f"{size_gb:.1f}G"
    elif size_bytes >= 1024 ** 2:
        size_mb = size_bytes / (1024 ** 2)
        return f"{size_mb:.0f}M" if size_mb >= 10 else f"{size_mb:.1f}M"
    elif size_bytes >= 1024:
        return f"{size_bytes / 1024:.1f}K"
    else:
        return f"{size_bytes}B"
