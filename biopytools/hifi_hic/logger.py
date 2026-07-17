"""
基因组组装日志管理模块|Genome Assembly Logger Module

遵循§2.3超算日志分离规范:INFO→stdout→.out,WARNING+→stderr→.err,DEBUG+→file
|Follows §2.3 job scheduler log separation: INFO→stdout→.out, WARNING+→stderr→.err, DEBUG+→file
"""

import logging
import sys
from pathlib import Path


class AssemblyLogger:
    """组装日志管理器|Assembly Logger Manager"""

    def __init__(self, output_dir, log_name: str = "assembly.log"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_dir: 日志目录|Log directory (自动创建|auto-created)
            log_name: 日志文件名|Log filename
        """
        self.output_dir = Path(output_dir)
        self.log_file = self.output_dir / log_name
        self._setup_logging()

    def _setup_logging(self):
        """设置日志(§2.3标准:stdout/stderr/file分离)|Setup logging (§2.3 standard)"""
        # 确保日志目录存在|Ensure log directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 独立logger(避免污染root,避免重复)|Independent logger (no root pollution, no duplicates)
        logger = logging.getLogger("hifi_hic")
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - INFO级别 → 超算捕获到.out文件
        # |stdout handler - INFO level → captured by scheduler to .out
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上 → 超算捕获到.err文件
        # |stderr handler - WARNING and above → captured by scheduler to .err
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # file handler - DEBUG及以上,追加模式(不删历史,支持断点续传排查)
        # |file handler - DEBUG and above, append mode (keep history for resume debugging)
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8', mode='a')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger
