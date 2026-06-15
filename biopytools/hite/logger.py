"""
HiTE 日志管理模块|HiTE Logging Management Module

提供统一的日志管理功能，支持HiTE和panHiTE
Provides unified logging management for both HiTE and panHiTE
"""

import logging
import sys
from pathlib import Path
from datetime import datetime


class HiteLoggerManager:
    """
    HiTE 日志管理器|HiTE Logger Manager

    为HiTE和panHiTE提供统一的日志记录功能
    Provides unified logging functionality for HiTE and panHiTE
    """

    def __init__(self, output_dir: str, log_prefix: str = "hite", log_level: str = "INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_dir: 输出目录路径|Output directory path
            log_prefix: 日志文件前缀|Log file prefix (default: "hite")
            log_level: 日志级别|Log level (default: "INFO")
        """
        self.output_dir = Path(output_dir)
        self.log_prefix = log_prefix
        self.log_level = log_level

        # 创建日志文件名|Create log filename
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.log_file = self.output_dir / f"{log_prefix}_{timestamp}.log"

        # 设置日志|Setup logging
        self.setup_logging()

    def setup_logging(self):
        """
        设置日志系统|Setup logging system

        配置日志格式和处理器
        Configure log format and handlers
        """
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        # 设置日志级别|Set log level
        level = getattr(logging, self.log_level.upper(), logging.INFO)

        # 创建处理器|Create handlers
        handlers = [logging.StreamHandler(sys.stdout)]
        handlers.append(logging.FileHandler(self.log_file))

        # 配置日志|Configure logging
        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers
        )

        self.logger = logging.getLogger(self.__class__.__name__)

        # 记录日志文件位置|Log log file location
        self.logger.info(f"日志文件|Log file: {self.log_file}")
        self.logger.info(f"日志级别|Log level: {self.log_level}")

    def get_logger(self):
        """
        获取日志器|Get logger

        Returns:
            logging.Logger: 配置好的日志器对象|Configured logger object
        """
        return self.logger
