"""reads2tree 工具函数|reads2tree utilities

fastq 目录扫描、双端配对、reads 合并(cat/BBMerge)、样本映射与 waster 输入清单。
复用 genome2tree 的映射/清单工具与 eviann 的双端配对逻辑,不造第二套。
|fastq dir scanning, paired-end grouping, read merging (cat/BBMerge), sample
map parsing and waster input list generation. Reuses genome2tree utilities and
eviann's pairing logic.

日志管理器(三 handler:stdout/stderr/file)|Module logger (3 handlers)
"""

import gzip
import logging
import os
import shutil
import sys
from dataclasses import dataclass
from typing import List, Optional, Tuple

from ..eviann.calculator import pair_short_reads
from ..genome2tree.utils import (parse_samples_map, write_input_tsv,
                                 write_samples_map_tsv)

# fastq 后缀(小写匹配,大小写不敏感)|FASTQ extensions (case-insensitive)
FASTQ_EXTS = {".fq", ".fastq"}


class ModuleLogger:
    """模块日志管理器(三 handler:stdout/stderr/file)|Module logger (3 handlers)"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.logger = logging.getLogger("reads2tree")
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


@dataclass
class SampleFile:
    """单个样本合并后的 reads 文件|One sample's merged reads file"""
    stem: str        # 样本名|sample name
    path: str        # 合并后的明文 fastq 绝对路径|merged plaintext fastq path


def split_fastq_filename(filename: str) -> Optional[Tuple[str, bool]]:
    """拆 fastq 文件名|Split a FASTQ filename.

    后缀匹配大小写不敏感;支持 .gz 双后缀;只认 fastq 系后缀。
    |Case-insensitive; .gz double suffix; FASTQ extensions only.

    Returns:
        (stem, is_gz);非 fastq 文件返回 None|(stem, is_gz); None otherwise.
    """
    lowered = filename.lower()
    is_gz = lowered.endswith(".gz")
    if is_gz:
        lowered = lowered[:-3]
    for ext in FASTQ_EXTS:
        if lowered.endswith(ext):
            cut = len(ext) + (3 if is_gz else 0)
            return filename[:-cut], is_gz
    return None


def scan_fastq_dir(
        input_dir: str
) -> Tuple[List[str], List[str], List[Tuple[str, List[str], List[str]]], List[str]]:
    """递归扫描输入目录并配对双端 reads|Recursively scan input dir and group paired reads

    支持多级子目录(如 raw/、clean/);非 fastq 文件与空目录忽略(容错不报错)。
    配对规则复用 eviann(支持 _R1/_R2、_1/_2、.R1.、read1/read2、_R1_001、
    _1_clean.fq.gz 等形态);同一样本多 lane 归入同组;缺一侧的归单端。
    |Multi-level subdirs supported (e.g. raw/, clean/); non-fastq files ignored.
    Pairing reuses eviann rules (incl. _1_clean.fq.gz); multi-lane samples group
    together; one-sided pairs fall back to single-end.

    Raises:
        ValueError: 同一样本出现在多个子目录(如 raw/ 与 clean/ 都有 SRRxxx,
        避免原始+clean 混拼)|same sample in multiple subdirs (avoid mixing raw+clean)

    Returns:
        (fastq_paths, ignored, paired, singles)
        paired: (样本名|sample, R1 列表|R1 list, R2 列表|R2 list)
    """
    fastq_paths: List[str] = []
    ignored: List[str] = []
    for root, _dirs, files in os.walk(input_dir):
        for name in sorted(files):
            path = os.path.join(root, name)
            if split_fastq_filename(name) is None:
                ignored.append(os.path.relpath(path, input_dir))
            else:
                fastq_paths.append(os.path.abspath(path))
    paired, singles = pair_short_reads(fastq_paths)

    # 冲突检测:同一样本跨多个子目录 → 报错|conflict: sample in multiple subdirs
    sample_dirs: dict = {}
    for sample, r1s, r2s in paired:
        sample_dirs.setdefault(sample, set()).update(
            _relative_dir(input_dir, p) for p in r1s + r2s)
    for p in singles:
        stem = split_fastq_filename(os.path.basename(p))
        sample = stem[0] if stem else os.path.basename(p)
        sample_dirs.setdefault(sample, set()).add(_relative_dir(input_dir, p))
    conflicts = {s: sorted(dirs) for s, dirs in sample_dirs.items() if len(dirs) > 1}
    if conflicts:
        detail = "; ".join(f"{s} ({', '.join(dirs)})"
                           for s, dirs in sorted(conflicts.items()))
        raise ValueError(
            f"以下样本出现在多个子目录|Samples found in multiple subdirs: {detail}。"
            f"请将输入目录指向具体数据目录(如 clean/)|point -i at the actual "
            f"data dir (e.g. clean/)")

    return fastq_paths, ignored, paired, singles


def _relative_dir(input_dir: str, path: str) -> str:
    """文件所在子目录(相对输入目录)|Subdir of a file relative to input dir"""
    rel = os.path.relpath(os.path.dirname(path), input_dir)
    return rel if rel != "." else "(根目录|root)"


def concat_reads(files: List[str], dest: str,
                 logger: Optional[logging.Logger] = None) -> bool:
    """流式拼接多个 fastq 到单个明文文件|Stream-concat fastq files into one plaintext file

    gz 文件逐流解压、明文直接复制;dest 已存在则跳过(断点续传)。
    |Gzipped files are stream-decompressed, plaintext copied; skip if dest exists.

    Returns:
        bool: 成功|success
    """
    if os.path.exists(dest):
        if logger:
            logger.info(f"跳过已完成步骤|Skipping: 已合并|already merged: {dest}")
        return True
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if logger:
        logger.info(f"合并 reads|Merging reads -> {dest} "
                    f"({len(files)} 文件|files)")
    try:
        with open(dest, "wb") as out:
            for f in files:
                opener = gzip.open if f.lower().endswith(".gz") else open
                with opener(f, "rb") as fin:
                    shutil.copyfileobj(fin, out)
        return True
    except OSError as e:
        if logger:
            logger.error(f"合并失败|Merge failed: {e}")
        return False


def merge_paired_reads(
        sample: str, r1_files: List[str], r2_files: List[str],
        dest_dir: str, bbmerge_path: str,
        logger: Optional[logging.Logger] = None) -> Optional[str]:
    """用 BBMerge 合并重叠双端 reads|Merge overlapping paired reads with BBMerge

    每对 lane 单独 BBMerge(-in1/-in2/-out),产物再拼接为 {sample}.fq。
    全部 -outu 未合并 reads 丢弃(WASTER 只吃合并代表)。
    |Each lane pair is merged separately; outputs are concatenated into
    {sample}.fq. Unmerged reads (-outu) are discarded.

    Returns:
        str: 合并产物路径|merged path, None 表示失败|None on failure
    """
    import subprocess

    dest = os.path.join(dest_dir, f"{sample}.fq")
    if os.path.exists(dest):
        if logger:
            logger.info(f"跳过已完成步骤|Skipping: 已合并|already merged: {dest}")
        return dest
    os.makedirs(dest_dir, exist_ok=True)

    parts = []
    for i, (r1, r2) in enumerate(zip(r1_files, r2_files)):
        part = os.path.join(dest_dir, f"{sample}_lane{i}_merged.fq")
        cmd = [bbmerge_path, "-in1", r1, "-in2", r2, "-out", part]
        if logger:
            logger.info(f"命令|Command: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, shell=False, check=False,
                                    capture_output=True, text=True)
        except FileNotFoundError:
            if logger:
                logger.error(f"BBMerge 未找到|bbmerge.sh not found: {bbmerge_path}")
            return None
        if result.returncode != 0:
            if logger:
                logger.error(f"BBMerge 失败|BBMerge failed (rc={result.returncode}): "
                             f"{(result.stderr or result.stdout or '').strip()[:500]}")
            return None
        if not os.path.exists(part):
            if logger:
                logger.warning(f"BBMerge 未产出|BBMerge produced nothing for {r1}")
            continue
        parts.append(part)

    if not parts:
        if logger:
            logger.error(f"样本无可合并 reads|No mergeable reads for {sample}")
        return None
    if not concat_reads(parts, dest, logger):
        return None
    return dest
