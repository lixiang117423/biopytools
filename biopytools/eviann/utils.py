"""
EviAnn工具函数模块|EviAnn Utilities Module
"""

import logging
import sys


class ModuleLogger:
    """模块日志管理器|Module Logger Manager

    使用独立命名空间避免污染root|Named logger to avoid root pollution
    """

    def __init__(self, log_file=None, log_level="INFO"):
        self.log_file = log_file
        self.logger = self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        logger = logging.getLogger('biopytools.eviann')
        logger.setLevel(level)
        logger.handlers.clear()
        logger.propagate = False

        formatter = logging.Formatter(log_format, date_format)

        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

        if self.log_file:
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        return logger
