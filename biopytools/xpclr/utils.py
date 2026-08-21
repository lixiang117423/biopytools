"""xpclr 工具函数|xpclr utilities.

日志管理器 + 大数字格式化 + xpclr 版本探测(CLI 无 --version,从包源码读)。
命令执行统一走 common/conda_runner.CommandRunner(§13.2,完整命令记 INFO)。
|Logger + number formatting + xpclr version probing; command execution via
the common CommandRunner layer.
"""
from __future__ import annotations

import glob
import logging
import os
import re
import sys
from pathlib import Path
from typing import Optional

LOG_FORMAT = "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s"
LOG_DATEFMT = "%Y-%m-%d %H:%M:%S"


class ModuleLogger:
    """模块日志管理器(三 handler: stdout/stderr/file)|Module logger (3 handlers)."""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.logger = logging.getLogger("xpclr")
        self.logger.handlers.clear()
        self.logger.propagate = False
        self.logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        fmt = logging.Formatter(LOG_FORMAT, LOG_DATEFMT)
        sh = logging.StreamHandler(sys.stdout)   # INFO+ → 超算 .out|stdout
        sh.setLevel(logging.INFO)
        sh.setFormatter(fmt)
        self.logger.addHandler(sh)
        eh = logging.StreamHandler(sys.stderr)   # WARNING+ → 超算 .err|stderr
        eh.setLevel(logging.WARNING)
        eh.setFormatter(fmt)
        self.logger.addHandler(eh)
        if log_file:                              # 全级别 → 文件|file
            Path(log_file).parent.mkdir(parents=True, exist_ok=True)
            fh = logging.FileHandler(log_file)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(fmt)
            self.logger.addHandler(fh)

    def get_logger(self) -> logging.Logger:
        """返回 logger|Return logger."""
        return self.logger


def format_number(num) -> str:
    """大于1百万用 M 单位保留2位小数(§5.3)|Format >1M numbers in M units."""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)


def detect_xpclr_version(xpclr_path: str) -> str:
    """读 xpclr 包 __init__.py 的 __version__|Read __version__ from package source.

    xpclr CLI 无 --version 参数;从 <env>/lib/python*/site-packages/xpclr/__init__.py
    正则抽取,失败返回 unknown(不抛错,版本探测不阻塞流程)。
    |No --version flag on the CLI; regex from site-packages, "unknown" on failure.
    """
    env_root = os.path.dirname(os.path.dirname(os.path.abspath(xpclr_path)))
    candidates = glob.glob(
        os.path.join(env_root, "lib", "python*", "site-packages", "xpclr", "__init__.py"))
    for cand in candidates:
        try:
            m = re.search(r'__version__\s*=\s*["\']([^"\']+)', Path(cand).read_text())
            if m:
                return m.group(1)
        except OSError:
            continue
    return "unknown"
