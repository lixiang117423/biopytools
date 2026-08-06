"""
FOF生成工具函数模块|FOF Generation Utility Functions Module
"""

import logging
import sys
from pathlib import Path
from typing import List, Optional

# 压缩后缀|Compression suffixes
_COMPRESSION_EXTS = ('.gz', '.bz2', '.xz', '.zst')
# 常见序列/比对数据扩展|Common sequence/alignment data extensions
_DATA_EXTS = (
    '.fastq', '.fq', '.fasta', '.fa', '.fna', '.faa', '.ffn',
    '.bam', '.sam', '.cram',
    '.vcf', '.bcf',
    '.bed', '.gtf', '.gff3', '.gff',
)


class FofLogger:
    """FOF生成日志管理器(三 handler:stdout/stderr/file)|FOF Logger (3 handlers)"""

    def __init__(self, output_dir: Path, log_name: str = "fof.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 删除已存在的日志文件|Remove existing log file
        if self.log_file.exists():
            self.log_file.unlink()

        # 标准日志格式|Standard log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 命名 logger + 清空旧 handler + 不向 root 传播(避免重复)|Named logger, clear handlers, no propagation
        self.logger = logging.getLogger("FOF")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False

        # stdout handler - INFO级别 → 超算 .out|stdout - INFO → HPC .out
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # stderr handler - WARNING及以上 → 超算 .err|stderr - WARNING+ → HPC .err
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        # 文件 handler - 所有级别 → 本地完整日志|file - all levels → local full log
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(stderr_handler)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def extract_sample_name(filepath: Path, suffixes: Optional[List[str]] = None) -> str:
    """
    从文件名提取样品名(剥离数据扩展,含压缩后缀)|Extract sample name (strip data ext incl. compression)

    优先级|Priority:
        1. 若提供 suffixes 且文件名以某个后缀结尾,精确剥离该后缀(最长匹配)
           |If suffixes given and filename ends with one, strip it (longest match)
        2. 否则:先剥压缩后缀(.gz/.bz2/...),再剥一层已知数据扩展(.fastq/.bam/...)
           |Else: strip compression suffix then one known data extension
        3. 兜底:剥最后一层扩展(保持旧行为)|Fallback: strip last extension (legacy)

    Args:
        filepath: 文件路径|File path
        suffixes: 已知过滤后缀列表(如 ['.fastq.gz', '.fq.gz'])|Known filter suffixes

    Returns:
        样品名|Sample name

    Examples:
        >>> extract_sample_name(Path('/data/sample1_R1.fastq.gz'))
        'sample1_R1'
        >>> extract_sample_name(Path('/data/sample1.bam'))
        'sample1'
        >>> extract_sample_name(Path('/data/sample1_R1.fq.gz'), suffixes=['.fq.gz'])
        'sample1_R1'
    """
    filename = filepath.name
    name_lower = filename.lower()

    # 1. 已知后缀精确剥离(最长匹配优先)|Exact strip of known suffix (longest first)
    if suffixes:
        for suf in sorted(suffixes, key=len, reverse=True):
            if suf and name_lower.endswith(suf.lower()):
                return filename[:-len(suf)] or filename

    # 2. 先剥压缩后缀,再剥已知数据扩展|Strip compression suffix then one known data ext
    for comp in _COMPRESSION_EXTS:
        if name_lower.endswith(comp):
            filename = filename[:-len(comp)]
            name_lower = filename.lower()
            break
    for ext in _DATA_EXTS:
        if name_lower.endswith(ext):
            return filename[:-len(ext)] or filename

    # 3. 兜底:剥最后一层扩展|Fallback: strip last extension
    dot_pos = filename.rfind('.')
    if dot_pos > 0:
        return filename[:dot_pos]
    return filename


def format_number(num: int) -> str:
    """格式化数字|Format number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)
