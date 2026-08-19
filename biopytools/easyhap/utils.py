"""easyhap 工具函数|easyhap utilities"""
import gzip
import logging
import re
import sys
from typing import List, Optional, Tuple

# 区域字符串格式 chr:start-end|region string format
REGION_RE = re.compile(r"^([^:]+):(\d+)-(\d+)$")


class ModuleLogger:
    """模块日志管理器(三 handler: stdout/stderr/file)|Module logger"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.logger = logging.getLogger("easyhap")
        self.logger.handlers.clear()
        self.logger.propagate = False
        self.logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        fmt = logging.Formatter(
            "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            "%Y-%m-%d %H:%M:%S")
        sh = logging.StreamHandler(sys.stdout)   # INFO+ → 超算 .out
        sh.setLevel(logging.INFO)
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


def parse_region_string(region: str) -> Tuple[str, int, int]:
    """解析 chr:start-end|Parse chr:start-end

    Returns:
        (chrom, start, end); 格式错误或 start>end 抛 ValueError
        |raises ValueError on bad format or start>end
    """
    m = REGION_RE.match(region.strip())
    if not m:
        raise ValueError(
            f"区域格式错误(需 chr:start-end)|Bad region format (need chr:start-end): {region!r}")
    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    if start > end:
        raise ValueError(f"区域起始大于终止|Region start > end: {region!r}")
    return chrom, start, end


def read_region_file(region_file: str) -> List[Tuple[str, int, int]]:
    """解析批量区域文件(TAB 三列 chr start end, 无表头)|Parse batch region file

    支持 .gz|gz supported; `#` 注释与空行忽略; 非法行 ValueError 带行号
    |comments/blank skipped; malformed lines raise with line number
    """
    regions: List[Tuple[str, int, int]] = []
    opener = gzip.open if region_file.endswith(".gz") else open
    with opener(region_file, "rt") as fh:
        for line_no, raw in enumerate(fh, 1):
            line = raw.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            if "\t" not in line:
                raise ValueError(
                    f"区域文件第 {line_no} 行须 TAB 分隔|region file line {line_no} "
                    f"must be TAB-delimited: {line!r}")
            parts = line.split("\t")
            if len(parts) < 3 or not all(p.strip() for p in parts[:3]):
                raise ValueError(
                    f"区域文件第 {line_no} 行需三列 chr<TAB>start<TAB>end|region file "
                    f"line {line_no} needs chr<TAB>start<TAB>end: {line!r}")
            try:
                start, end = int(parts[1].strip()), int(parts[2].strip())
            except ValueError:
                raise ValueError(
                    f"区域文件第 {line_no} 行起止非整数|region file line {line_no} "
                    f"start/end not integer: {line!r}")
            if start > end:
                raise ValueError(
                    f"区域文件第 {line_no} 行起始大于终止|region file line {line_no} "
                    f"start > end: {line!r}")
            regions.append((parts[0].strip(), start, end))
    return regions


def region_label(chrom: str, start: int, end: int) -> str:
    """输出前缀标签(与上游一致)|Output prefix label (same as upstream)

    上游: sanitize_filename(f"{chrom}_{start}_{end}"), 非 [A-Za-z0-9_.-] 字符替换为 _
    """
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", f"{chrom}_{start}_{end}")
