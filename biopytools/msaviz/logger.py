"""
日志管理模块|Logger Management Module
"""

import logging
import sys
from pathlib import Path


class MsaVizLogger:
    """MSA可视化日志管理器|MSA Visualization Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Parameters
        ----------
        log_file : str or Path, optional
            日志文件路径|Log file path
        log_level : str, optional
            日志级别|Log level (default: "INFO")
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """
        设置日志|Setup logging

        Parameters
        ----------
        log_level : str
            日志级别|Log level
        """
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """
        获取日志器|Get logger

        Returns
        -------
        logger : logging.Logger
            日志器对象|Logger object
        """
        return self.logger
