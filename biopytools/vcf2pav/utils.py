"""vcf2pav 工具函数|vcf2pav utility functions

日志管理器、GT解码、SV ID生成。
|Logger, GT decoding, SV ID generation.
"""
import logging
import sys
from typing import Optional


class ModuleLogger:
    """模块日志管理器(三 handler)|Module logger (3 handlers)"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.logger = logging.getLogger("vcf2pav")
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


def decode_gt(gt_field: str) -> Optional[int]:
    """解码 GT 字段为 0/1|Decode GT field to 0/1.

    Args:
        gt_field: GT 字段值 (如 "1/1", "./.", "0/1", "1|0")
                  |GT field value (e.g. "1/1", "./.", "0/1", "1|0")

    Returns:
        1=存在(present), 0=缺失(absent), None=无法解析(unparseable)
    """
    if not gt_field or gt_field == ".":
        return 0
    if gt_field in ("./.", ".|."):
        return 0
    if gt_field in ("0/0", "0|0"):
        return 0
    if gt_field in ("1/1", "0/1", "1/0", "1|1", "0|1", "1|0"):
        return 1
    # 尝试通用解析: 等位基因须为数字/点号,含非数字视为无法解析
    # |Generic parse: alleles must be numeric/dot, non-numeric → unparseable
    alleles = gt_field.replace("/", "|").split("|")
    if not all(a in ("0", ".", "") or a.isdigit() for a in alleles):
        return None
    if any(a not in ("0", ".", "") for a in alleles):
        return 1
    # 如果全是 0 或 . → absent
    return 0


def make_sv_id(chrom: str, start: int, end: int, svtype: str,
               chr2: str = "") -> str:
    """生成 SV 唯一标识符|Generate unique SV identifier.

    Format:
        BND/TRA 跨染色体: {CHROM}_{START}-{CHR2}_{END}_{SVTYPE}
        其他: {CHROM}_{START}-{END}_{SVTYPE}
    """
    if svtype in ("BND", "TRA") and chr2 and chr2 != chrom:
        return f"{chrom}_{start}-{chr2}_{end}_{svtype}"
    return f"{chrom}_{start}-{end}_{svtype}"
