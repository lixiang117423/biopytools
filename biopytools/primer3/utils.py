"""
Primer3引物设计工具函数模块|Primer3 Primer Design Utility Functions Module
"""

import logging
import re
import sys
from pathlib import Path

# conda 调用统一走公共层(禁本地副本, §13.3)
# |Conda invocation from the common layer (no local copies, §13.3)
from ..common.conda_runner import CommandRunner, build_conda_command, get_conda_env

__all__ = [
    'Primer3Logger',
    'CommandRunner',
    'build_conda_command',
    'get_conda_env',
    'format_sequence_id',
    'format_number',
]


class Primer3Logger:
    """Primer3引物设计日志管理器(named logger, §2.3 三通道分离)
    |Primer3 Logger Manager (named logger, §2.3 three-stream separation)

    stdout(<=INFO)→.out + stderr(>=WARNING)→.err + 文件(DEBUG 全量);
    named logger 不污染 root|named logger, no root pollution
    """

    def __init__(self, log_dir: Path, log_name: str = "primer3_design.log"):
        self.log_dir = Path(log_dir)
        self.log_file = self.log_dir / log_name
        self.logger = self._setup_logging()

    def _setup_logging(self) -> logging.Logger:
        """配置 named logger 与三 handler|Configure named logger with 3 handlers."""
        self.log_dir.mkdir(parents=True, exist_ok=True)
        # 重跑时重建日志文件, 避免新旧运行混写|Fresh log per run
        if self.log_file.exists():
            self.log_file.unlink()

        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        logger = logging.getLogger('biopytools.primer3')
        logger.handlers.clear()
        logger.propagate = False
        logger.setLevel(logging.DEBUG)

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        # INFO 及以下才进 stdout, WARNING 及以上只走 stderr(§2.3)
        # |Only <=INFO goes to stdout; >=WARNING goes to stderr only (§2.3)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        return logger

    def get_logger(self) -> logging.Logger:
        """获取日志器|Get logger"""
        return self.logger


def format_sequence_id(seq_id: str) -> str:
    """
    格式化序列ID，确保可以作为Primer3输入|Format sequence ID for Primer3 input

    Primer3要求SEQUENCE_ID不能包含特殊字符和空格|Primer3 requires SEQUENCE_ID without special chars and spaces

    Args:
        seq_id: 原始序列ID|Original sequence ID

    Returns:
        格式化后的序列ID|Formatted sequence ID
    """
    formatted = re.sub(r'[^\w\-.]', '_', seq_id)
    return formatted


def format_number(num: int) -> str:
    """大数字格式化: ≥1百万用 M 单位保留2位(§5.3)|Large numbers use M unit (§5.3)

    Args:
        num: 数字|Number

    Returns:
        格式化字符串|Formatted string
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)
