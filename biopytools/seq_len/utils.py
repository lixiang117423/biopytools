"""seq_len 工具函数|seq_len utilities (logger, streaming FASTA length, N50, summary)"""

import gzip
import logging
import sys
from typing import Dict, Iterator, List, Optional, Tuple


class SeqLenLogger:
    """三分日志:stdout=INFO / stderr=WARNING+ / file=DEBUG+|3-way split logger"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = 'INFO',
                 verbose: bool = False):
        self.log_file = log_file
        self.log_level = logging.DEBUG if verbose else getattr(
            logging, str(log_level).upper(), logging.INFO)

    def get_logger(self) -> logging.Logger:
        """构建并返回配置好的 named logger|Build and return a configured logger"""
        logger = logging.getLogger('seq_len')
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False  # 避免向 root 重复输出|no duplicate output to root

        fmt = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')

        stdout_h = logging.StreamHandler(sys.stdout)
        stdout_h.setLevel(self.log_level)
        stdout_h.setFormatter(fmt)
        logger.addHandler(stdout_h)

        stderr_h = logging.StreamHandler(sys.stderr)
        stderr_h.setLevel(logging.WARNING)
        stderr_h.setFormatter(fmt)
        logger.addHandler(stderr_h)

        if self.log_file:
            file_h = logging.FileHandler(self.log_file)
            file_h.setLevel(logging.DEBUG)
            file_h.setFormatter(fmt)
            logger.addHandler(file_h)
        return logger


def format_number(num: int) -> str:
    """大数字(≥1M)用 M 单位保留2位小数|Big numbers (≥1M) in M unit, 2 decimals"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)


def _open_text(path: str):
    """透明打开 .gz 文本|Transparently open .gz text"""
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def iter_seq_lengths(path: str) -> Iterator[Tuple[str, int]]:
    """流式产出 (seq_id, length),不把整条序列载入内存|Stream (seq_id, length) without loading whole sequences

    seq_id 取 header 去掉 '>' 后的首个空白 token;空序列记长度 0。
    |seq_id = first whitespace token of header; empty sequences count as 0.
    """
    seq_id: Optional[str] = None
    length = 0
    with _open_text(path) as f:
        for line in f:
            line = line.rstrip('\n').rstrip('\r')
            if line.startswith('>'):
                if seq_id is not None:
                    yield (seq_id, length)
                header = line[1:].strip()
                seq_id = header.split()[0] if header else ''
                length = 0
            elif line:
                length += len(line.strip())
        if seq_id is not None:
            yield (seq_id, length)


def compute_n50(lengths: List[int]) -> int:
    """计算 N50:长度降序累加,达到总长一半时的长度值|N50: length at which cumulative ≥ half of total"""
    if not lengths:
        return 0
    half = sum(lengths) / 2
    cum = 0
    for length in sorted(lengths, reverse=True):
        cum += length
        if cum >= half:
            return length
    return 0  # 不可达,防御性返回|unreachable, defensive


def compute_summary(lengths: List[int]) -> Dict[str, int]:
    """汇总统计:序列数/总长/N50/最短/最长/平均|Summary: count/total/N50/min/max/mean"""
    if not lengths:
        return {"num_seqs": 0, "total_length": 0, "n50": 0,
                "min_length": 0, "max_length": 0, "mean_length": 0}
    total = sum(lengths)
    num = len(lengths)
    return {
        "num_seqs": num,
        "total_length": total,
        "n50": compute_n50(lengths),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": round(total / num, 2),
    }
