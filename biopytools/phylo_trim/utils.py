"""phylo_trim 工具函数模块|phylo_trim Utility Functions Module"""

import logging
import sys


class PhyloTrimLogger:
    """phylo_trim 日志管理器|phylo_trim Logger Manager

    用 named logger "phylo_trim" + propagate=False,与 mafft_fasttree 的 root logger 隔离
    |Named logger "phylo_trim" with propagate=False, isolated from mafft_fasttree's root logger
    """

    def __init__(self, log_file=None, verbose=False):
        self.log_file = log_file
        self.setup_logging(verbose)

    def setup_logging(self, verbose):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        formatter = logging.Formatter(log_format, datefmt=date_format)

        level = logging.DEBUG if verbose else logging.INFO

        logger = logging.getLogger("phylo_trim")
        logger.setLevel(level)
        logger.handlers.clear()
        logger.propagate = False

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        if self.log_file:
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger
