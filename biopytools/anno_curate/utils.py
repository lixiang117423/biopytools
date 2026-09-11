"""anno_curate工具函数|anno_curate Utilities"""

import logging
import os
import sys


class AnnoCurateLogger:
    """anno_curate日志管理器|anno_curate Logger Manager

    stdout(<=INFO) + stderr(>=WARNING) + 三文件:
    anno_curate.log(全量) anno_curate.out.log(<=INFO) anno_curate.err.log(>=WARNING)
    """

    def __init__(self, logs_dir: str, log_level: str = 'INFO'):
        self.log_file = os.path.join(logs_dir, 'anno_curate.log')
        self.out_log_file = os.path.join(logs_dir, 'anno_curate.out.log')
        self.err_log_file = os.path.join(logs_dir, 'anno_curate.err.log')
        os.makedirs(logs_dir, exist_ok=True)
        self.logger = self._setup_logging(log_level)

    def _setup_logging(self, log_level: str) -> logging.Logger:
        """设置日志(named logger)|Setup logging (named logger)"""
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')
        level = getattr(logging, log_level.upper(), logging.INFO)

        logger = logging.getLogger('biopytools.anno_curate')
        logger.handlers.clear()
        logger.propagate = False
        logger.setLevel(logging.DEBUG)

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(level)
        stdout_handler.addFilter(lambda r: r.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        specs = [(self.log_file, None),
                 (self.out_log_file, lambda r: r.levelno <= logging.INFO),
                 (self.err_log_file, lambda r: r.levelno >= logging.WARNING)]
        for path, level_filter in specs:
            handler = logging.FileHandler(path, encoding='utf-8')
            handler.setLevel(logging.DEBUG)
            if level_filter:
                handler.addFilter(level_filter)
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        return logger

    def get_logger(self) -> logging.Logger:
        """获取日志器|Get logger"""
        return self.logger


def format_number(num: int) -> str:
    """大于1百万用M单位2位小数(§5.3)|Numbers >1M use M unit, 2 decimals"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)
