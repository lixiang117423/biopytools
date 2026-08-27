"""phy2fa 工具函数|phy2fa utilities

日志管理器、Phylip 解析与 FASTA 写出。
|Logger manager, Phylip parsing and FASTA writing.
"""

import gzip
import logging
import os
import sys
from typing import Optional


class ModuleLogger:
    """模块日志管理器(三 handler:stdout/stderr/file)|Module logger (3 handlers)"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.logger = logging.getLogger("phy2fa")
        self.logger.handlers.clear()
        self.logger.propagate = False
        self.logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
        fmt = logging.Formatter(
            "%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            "%Y-%m-%d %H:%M:%S")
        sh = logging.StreamHandler(sys.stdout)
        sh.setLevel(logging.INFO)
        sh.addFilter(lambda r: r.levelno <= logging.INFO)
        sh.setFormatter(fmt)
        self.logger.addHandler(sh)
        eh = logging.StreamHandler(sys.stderr)
        eh.setLevel(logging.WARNING)
        eh.setFormatter(fmt)
        self.logger.addHandler(eh)
        if log_file:
            fh = logging.FileHandler(log_file)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(fmt)
            self.logger.addHandler(fh)

    def get_logger(self) -> logging.Logger:
        return self.logger


def open_text(path: str):
    """按后缀透明打开 gz/明文文本|Transparent open for gz or plain text"""
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_header(line: str) -> tuple:
    """解析 Phylip 头行 (n_tax n_char)|Parse the Phylip header line

    Raises:
        ValueError: 头行格式非法|malformed header line
    """
    parts = line.split()
    if len(parts) != 2:
        raise ValueError(
            f"Phylip 头行格式错误(需 样本数 序列长度)|malformed header "
            f"(need n_tax n_char): {line!r}")
    try:
        n_tax, n_char = int(parts[0]), int(parts[1])
    except ValueError:
        raise ValueError(f"Phylip 头行含非整数|non-integer header: {line!r}")
    if n_tax <= 0 or n_char <= 0:
        raise ValueError(
            f"样本数/序列长度必须为正|n_tax/n_char must be positive: {line!r}")
    return n_tax, n_char


def parse_phy(lines, logger=None) -> dict:
    """解析 Phylip 数据行(已跳过首行头)→ {样本名: 完整序列}
    |Parse Phylip data rows (header already consumed) into records

    规则极简:每行按空白拆分,第一个字段为样本名,其余拼接为序列;
    同名再现追加片段(interleaved 续块)。
    |Dead simple: split each row on whitespace; first field = sample name,
    rest joined = sequence; same-name rows append (interleaved chunks).

    Args:
        lines: 迭代器/列表,已消费头行的剩余行|iterable after the header

    Returns:
        dict: 样本名→完整序列|sample name → full sequence
    """
    records = {}
    for raw in lines:
        s = raw.strip()
        if not s:
            continue
        parts = s.split()
        name = parts[0]
        seq = "".join(parts[1:]).upper()
        records[name] = records.get(name, "") + seq
    return records


def write_fasta(records: dict, out_path: str, line_width: int = 60,
                logger=None) -> None:
    """写 FASTA(换行宽度可选)|Write FASTA with optional wrapping"""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            if line_width and line_width > 0:
                for i in range(0, len(seq), line_width):
                    fh.write(seq[i:i + line_width] + "\n")
            else:
                fh.write(seq + "\n")
    if logger:
        logger.info(f"FASTA 已写入|FASTA written: {out_path} "
                    f"({len(records)} 条|records)")
