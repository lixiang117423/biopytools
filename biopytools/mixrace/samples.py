"""mixrace 样本配对识别|mixrace sample-pair discovery.

扫描 fastq 目录,按后缀配对 R1/R2。支持 _1/_2 与 _R1/_R2 两套命名(含 .clean. 变体)。
边界正确:用"匹配后缀后剥离"而非 glob 通配,S1_1.fq.gz 不会吞掉 S10_1.fq.gz。
|Scan a fastq dir and pair R1/R2 by suffix. Supports both _1/_2 and _R1/_R2 naming
(with .clean. variants). Boundary-correct: strips the matched suffix (not glob-wildcard),
so S1_1.fq.gz does not swallow S10_1.fq.gz.
"""
import os
import re
from typing import Dict, List

# (R1 后缀, R2 后缀),长在前避免误匹配|(R1 suffix, R2 suffix); longest first
_R1_R2_PAIRS = [
    ("_1.clean.fq.gz", "_2.clean.fq.gz"),
    ("_R1.clean.fq.gz", "_R2.clean.fq.gz"),
    ("_1.fastq.gz", "_2.fastq.gz"),
    ("_R1.fastq.gz", "_R2.fastq.gz"),
    ("_1.fq.gz", "_2.fq.gz"),
    ("_R1.fq.gz", "_R2.fq.gz"),
]


def _natural_key(name: str):
    """自然排序键:数字段按数值比较(Pb2 < Pb10,而非字典序 Pb10 < Pb2)。
    |Natural-sort key: digit runs compared numerically (Pb2 < Pb10)."""
    return [int(tok) if tok.isdigit() else tok
            for tok in re.split(r"(\d+)", name)]


def discover_samples(fastq_dir: str) -> List[Dict[str, str]]:
    """发现成对样本|Discover paired samples.

    Args:
        fastq_dir: 含 R1/R2 fastq.gz 的目录|dir with R1/R2 fastq.gz

    Returns:
        [{"sample": str, "r1": abs_path, "r2": abs_path}, ...] 按样本名排序
        |list of dicts (sorted by sample name)

    Raises:
        FileNotFoundError: 某 R1 找不到对应 R2|R1 has no matching R2
    """
    samples: Dict[str, Dict[str, str]] = {}
    for fname in sorted(os.listdir(fastq_dir)):
        for r1_suffix, r2_suffix in _R1_R2_PAIRS:
            if fname.endswith(r1_suffix):
                sample = fname[:-len(r1_suffix)]
                if sample in samples:
                    break  # 已由更长后缀命中(如 _1.clean.fq.gz)|already matched by longer suffix
                r1 = os.path.realpath(os.path.join(fastq_dir, fname))
                r2 = os.path.realpath(os.path.join(fastq_dir, sample + r2_suffix))
                if not os.path.exists(r2):
                    raise FileNotFoundError(
                        f"样本|sample {sample}: R2 not found (expected {sample + r2_suffix})")
                samples[sample] = {"sample": sample, "r1": r1, "r2": r2}
                break  # 一个文件只匹配一个后缀|one file matches one suffix
    return [samples[k] for k in sorted(samples, key=_natural_key)]
