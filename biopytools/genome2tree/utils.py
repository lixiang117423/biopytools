"""genome2tree 工具函数|genome2tree utilities

日志管理器、输入目录扫描、gz 解压、样本映射解析、waster 输入文件生成。
|Logger, input dir scanning, gz decompression, sample map parsing, waster input generation.
"""
import logging
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


class ModuleLogger:
    """模块日志管理器(三 handler:stdout/stderr/file)|Module logger (3 handlers)"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.logger = logging.getLogger("genome2tree")
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


# 序列文件后缀(小写匹配,大小写不敏感)|Sequence extensions (lowercase match)
FASTA_EXTS = {".fa", ".fasta", ".fna"}
FASTQ_EXTS = {".fq", ".fastq"}


@dataclass
class SampleFile:
    """单个样本输入文件|One sample input file"""
    stem: str        # 样本名(文件名去后缀去.gz)|sample name (filename minus ext minus .gz)
    path: str        # 实际路径(解压后指向解压副本)|actual path (decompressed copy if .gz)
    seq_type: str    # "fasta"|"fastq"
    is_gz: bool      # 原始文件是否 gz 压缩|whether original file is gzipped


def split_sequence_filename(filename: str) -> Optional[Tuple[str, str, bool]]:
    """拆序列文件名|Split a sequence filename.

    后缀匹配大小写不敏感;支持 .gz 双后缀。
    |Case-insensitive extension match; supports .gz double suffix.

    Returns:
        (stem, seq_type, is_gz);非序列文件返回 None。
        |(stem, seq_type, is_gz); None for non-sequence files.
    """
    lowered = filename.lower()
    is_gz = lowered.endswith(".gz")
    if is_gz:
        lowered = lowered[:-3]
    for exts, seq_type in ((FASTA_EXTS, "fasta"), (FASTQ_EXTS, "fastq")):
        for ext in sorted(exts):
            if lowered.endswith(ext):
                cut = len(ext) + (3 if is_gz else 0)
                return filename[:-cut], seq_type, is_gz
    return None


def scan_input_dir(input_dir: str) -> Tuple[List[SampleFile], List[str]]:
    """扫描输入目录|Scan input directory.

    只认单层文件;子目录与非序列文件归入 ignored(容错不报错,偏好能跑多少跑多少)。
    |Flat files only; subdirs and non-sequence files go to ignored (tolerant).

    Returns:
        (samples, ignored_files) 均按文件名排序|(samples, ignored), both sorted.
    """
    samples: List[SampleFile] = []
    ignored: List[str] = []
    for name in sorted(os.listdir(input_dir)):
        path = os.path.join(input_dir, name)
        if not os.path.isfile(path):
            ignored.append(name)
            continue
        parsed = split_sequence_filename(name)
        if parsed is None:
            ignored.append(name)
        else:
            stem, seq_type, is_gz = parsed
            samples.append(SampleFile(stem=stem, path=os.path.abspath(path),
                                      seq_type=seq_type, is_gz=is_gz))
    return samples, ignored


def decompress_gz_samples(samples: List[SampleFile], uncompressed_dir: str,
                          logger: Optional[logging.Logger] = None) -> List[SampleFile]:
    """将 gz 样本解压到 uncompressed_dir|Decompress gz samples into uncompressed_dir.

    waster 要求非压缩输入;非 gz 样本原样返回;已有解压副本复用(断点续传语义)。
    |waster requires unzipped input; non-gz pass through; existing copies reused.
    """
    import gzip
    import shutil
    from dataclasses import replace

    os.makedirs(uncompressed_dir, exist_ok=True)
    out: List[SampleFile] = []
    for s in samples:
        if not s.is_gz:
            out.append(s)
            continue
        inner_ext = ".fa" if s.seq_type == "fasta" else ".fq"
        dst = os.path.join(uncompressed_dir, s.stem + inner_ext)
        if not os.path.exists(dst):
            if logger:
                logger.info(f"解压|Decompressing: {s.path} -> {dst}")
            with gzip.open(s.path, "rb") as fin, open(dst, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        out.append(replace(s, path=os.path.abspath(dst)))
    return out


def parse_samples_map(map_path: str) -> Dict[str, str]:
    """解析个体→物种映射|Parse individual-to-species map.

    每行 个体stem<TAB>物种名;`#` 注释与空行忽略;缺 TAB/空字段/重复个体名抛 ValueError。
    |Each line stem<TAB>species; `#` comments/blanks ignored; missing TAB,
    empty field, or duplicated key raises ValueError.
    """
    mapping: Dict[str, str] = {}
    with open(map_path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if "\t" not in line:
                raise ValueError(
                    f"映射第 {lineno} 行格式错误(需 个体<TAB>物种)|"
                    f"malformed map line {lineno} (need individual<TAB>species): {line!r}")
            ind, species = (p.strip() for p in line.split("\t", 1))
            if not ind or not species:
                raise ValueError(
                    f"映射第 {lineno} 行存在空字段|empty field in map line {lineno}: {line!r}")
            if ind in mapping:
                raise ValueError(f"映射中个体名重复|duplicated individual in map: {ind}")
            mapping[ind] = species
    return mapping


def write_input_tsv(samples: List[SampleFile], out_path: str) -> None:
    """写 waster 输入清单(个体名<TAB>绝对路径)|Write waster input list."""
    with open(out_path, "w") as fh:
        for s in samples:
            fh.write(f"{s.stem}\t{s.path}\n")


def write_samples_map_tsv(species_of: Dict[str, str], out_path: str) -> None:
    """写 waster -a 映射文件(个体→物种,按个体名排序)|Write waster -a mapping file."""
    with open(out_path, "w") as fh:
        for ind in sorted(species_of):
            fh.write(f"{ind}\t{species_of[ind]}\n")


def format_number(num: float) -> str:
    """大数字格式化(≥1M 用 M 单位保留2位小数)|Format numbers ≥1M in M unit."""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(int(num))
