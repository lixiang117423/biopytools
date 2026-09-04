"""vcf2deepbsa 工具函数|vcf2deepbsa utilities

日志管理器、VCF 逐行解析、大数字格式化。
|Logger, line-by-line VCF parsing, large-number formatting.
"""
import gzip
import logging
import sys
from typing import List, Optional


class ModuleLogger:
    """模块日志管理器(三 handler)|Module logger (3 handlers)"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.logger = logging.getLogger("vcf2deepbsa")
        self.logger.handlers.clear()
        self.logger.propagate = False
        self.logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        fmt = logging.Formatter(
            "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            "%Y-%m-%d %H:%M:%S")
        # stdout: INFO+ → .out
        sh = logging.StreamHandler(sys.stdout)
        sh.setLevel(logging.INFO)
        sh.setFormatter(fmt)
        self.logger.addHandler(sh)
        # stderr: WARNING+ → .err
        eh = logging.StreamHandler(sys.stderr)
        eh.setLevel(logging.WARNING)
        eh.setFormatter(fmt)
        self.logger.addHandler(eh)
        # file: all levels → file
        if log_file:
            fh = logging.FileHandler(log_file)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(fmt)
            self.logger.addHandler(fh)

    def get_logger(self) -> logging.Logger:
        """返回 logger|Return logger"""
        return self.logger


class VcfRecord:
    """单行 VCF 数据记录|Single VCF data record"""

    def __init__(self, line: str):
        info = line.split("\t")
        self.CHROM = info[0]
        self.POS = info[1]
        self.ID = info[2]
        self.REF = info[3]
        self.ALT = info[4]
        self.QUAL = info[5]
        self.FILTER = info[6]
        self.INFO = info[7]
        self.FORMAT = info[8].split(":")
        self.GT: List[dict] = []
        for sample_field in info[9:]:
            values = sample_field.split(":")
            self.GT.append(dict(zip(self.FORMAT, values)))


class VcfReader:
    """VCF 流式读取器(支持 gz)|Streaming VCF reader (.gz supported)"""

    def __init__(self, vcf_path: str):
        self.samples: List[str] = []
        if vcf_path.endswith(".gz"):
            self._fh = gzip.open(vcf_path, "rt")
        else:
            self._fh = open(vcf_path, "r")
        # 读到第一条数据行为止,收集样本名|Read up to first data line, collect samples
        self._pending: Optional[str] = None
        line = self._next_line()
        while line is not None and line.startswith("#"):
            if line.startswith("#CHROM"):
                self.samples = line.split("\t")[9:]
            line = self._next_line()
        self._pending = line

    def _next_line(self) -> Optional[str]:
        """读下一条非空行|Read next non-blank line"""
        for raw in self._fh:
            line = raw.rstrip("\n")
            if line:
                return line
        return None

    def __enter__(self) -> "VcfReader":
        return self

    def __exit__(self, *exc):
        self.close()

    def __iter__(self) -> "VcfReader":
        return self

    def __next__(self) -> VcfRecord:
        if self._pending is None:
            raise StopIteration
        line = self._pending
        self._pending = self._next_line()
        return VcfRecord(line)

    def close(self):
        """关闭文件句柄|Close file handle"""
        self._fh.close()


def format_number(num: int) -> str:
    """大数字格式化(>1M 用 M 单位)|Format large numbers (>1M as M unit)"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)
