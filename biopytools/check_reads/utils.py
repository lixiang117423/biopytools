"""check_reads 工具函数|check_reads utilities

日志管理器(三 handler:stdout/stderr/file)与 fastq 文件名识别。
|Logger manager (3 handlers) and fastq filename detection.
"""

import logging
import os
import sys
from typing import Optional, Tuple

from ..reads2tree.utils import split_fastq_filename


class ModuleLogger:
    """模块日志管理器(三 handler:stdout/stderr/file)|Module logger (3 handlers)"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.logger = logging.getLogger("check_reads")
        self.logger.handlers.clear()
        self.logger.propagate = False
        self.logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        fmt = logging.Formatter(
            "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            "%Y-%m-%d %H:%M:%S")
        sh = logging.StreamHandler(sys.stdout)   # INFO+ → 超算 .out
        sh.setLevel(logging.INFO)
        sh.addFilter(lambda r: r.levelno <= logging.INFO)
        sh.setFormatter(fmt)
        self.logger.addHandler(sh)
        eh = logging.StreamHandler(sys.stderr)   # WARNING+ → 超算 .err
        eh.setLevel(logging.WARNING)
        eh.setFormatter(fmt)
        self.logger.addHandler(eh)
        if log_file:                              # 全级别 → 文件
            fh = logging.FileHandler(log_file)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(fmt)
            self.logger.addHandler(fh)

    def get_logger(self) -> logging.Logger:
        """返回 logger|Return logger"""
        return self.logger


def collect_fastq_files(input_dirs, logger=None) -> tuple:
    """递归收集多个目录下的 fastq 文件|Recursively collect fastq files from dirs

    Returns:
        (fastq_paths, ignored): 文件路径列表与忽略清单(相对各输入目录)
        |(fastq paths, ignored relpaths)
    """
    from ..reads2tree.utils import split_fastq_filename as _sf

    fastq_paths = []
    ignored = []
    for d in input_dirs:
        for root, _dirs, files in os.walk(d):
            for name in sorted(files):
                path = os.path.join(root, name)
                if _sf(name) is None:
                    ignored.append(os.path.relpath(path, d))
                else:
                    fastq_paths.append(os.path.abspath(path))
    return sorted(fastq_paths), ignored


def detect_read_pair(name: str) -> Optional[Tuple[str, str]]:
    """从文件名检测双端配对(复用 eviann 规则)|Detect paired-end direction (eviann rules)

    Returns:
        (样本名, '1'|'2') 或 None(单端/无法配对)|(sample, mate) or None
    """
    from ..eviann.calculator import extract_pair_prefix
    return extract_pair_prefix(name)
