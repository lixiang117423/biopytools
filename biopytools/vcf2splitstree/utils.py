"""vcf2splitstree 工具函数|vcf2splitstree utilities

日志管理器、VCF→p-distance 矩阵转换(numpy 向量化)、SplitsTree6 距离 CSV 写出。
|Logger manager, VCF→p-distance conversion (numpy vectorized), SplitsTree6 CSV output.
"""

import gzip
import logging
import os
import sys
from typing import Optional

import numpy as np


class ModuleLogger:
    """模块日志管理器(三 handler:stdout/stderr/file)|Module logger (3 handlers)"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.logger = logging.getLogger("vcf2splitstree")
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


def read_vcf_gt_matrix(vcf_path: str):
    """从 VCF 读取样本基因型并计算 p-distance 矩阵(numpy 向量化)
    |Read VCF genotypes and compute p-distance matrix (numpy vectorized)

    仅处理双等位位点;缺失基因型在配对比较中跳过。
    向量化:每位点按联合基因型 3*a1+a2 编码(取值 0..8,-1=缺失),
    差异位点数与共享位点数均用矩阵乘法求得。
    |Biallelic sites only; missing skipped per-pair. Each site is encoded as
    combined genotype 3*a1+a2; diff/shared counts via matrix multiplication.

    Returns:
        (labels, matrix): 样本名列表与 n×n 距离矩阵|labels and n×n distance matrix
    """
    labels = []
    gt_rows = []
    miss_rows = []

    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                labels = parts[9:]
                continue
            parts = line.rstrip("\n").rstrip().split("\t")
            if len(parts) < 10:
                continue
            n = len(labels)
            gt = np.full(n, -1, dtype=np.int16)   # -1 = missing
            miss = np.zeros(n, dtype=bool)
            for k, field in enumerate(parts[9:]):
                g = field.split(":")[0]
                sep = "/" if "/" in g else ("|" if "|" in g else None)
                if sep is None:
                    miss[k] = True
                    continue
                x, y = g.split(sep)
                if x == "." or y == ".":
                    miss[k] = True
                    continue
                try:
                    gt[k] = 3 * int(x) + int(y)
                except ValueError:
                    miss[k] = True
            gt_rows.append(gt)
            miss_rows.append(miss)

    n = len(labels)
    if n == 0:
        raise ValueError("VCF 中未找到样本|No samples found in VCF header")
    if not gt_rows:
        raise ValueError("VCF 中未找到变异位点|No variant sites found in VCF")

    G = np.vstack(gt_rows).astype(np.int32)     # m × n, -1 = missing
    miss = np.vstack(miss_rows)                  # m × n
    ok = (~miss)

    # 共享位点数|shared count per pair (matrix multiplication)
    shared_counts = ok.astype(np.int64).T @ ok.astype(np.int64)

    # 差异位点数:按联合基因型值 0..8 分解,纯矩阵乘法
    diff_counts = np.zeros((n, n), dtype=np.int64)
    for v in range(9):
        v_cnt = (ok & (G == v)).astype(np.int64)
        other_cnt = (ok & (G != v)).astype(np.int64)
        diff_counts += v_cnt.T @ other_cnt       # i 为 v 且 j 非 v

    shared_f = shared_counts.astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        dist = np.where(shared_f > 0, diff_counts / shared_f, np.nan)
    dist = np.nan_to_num(dist, nan=0.0)
    np.fill_diagonal(dist, 0.0)
    return labels, dist


def write_distance_csv(labels, matrix, out_path: str) -> None:
    """写 SplitsTree6 兼容的距离 CSV|Write SplitsTree6-compatible distance CSV

    不带首行样本数:SplitsTree6 的 CSVReader 要求首行含 ≥4 个逗号(自动探测);
    每行 label + n 个数值。GUI 打开该 CSV 会自动识别为距离矩阵并跑 NeighborNet。
    |No leading count line (SplitsTree6 auto-detects by first-line commas);
    one line per sample: label + n values. SplitsTree6 GUI opens it as a
    distance matrix and auto-runs NeighborNet.
    """
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    n = len(labels)
    with open(out_path, "w") as fh:
        for i in range(n):
            row = ",".join(f"{v:.6f}" for v in matrix[i])
            fh.write(f"{labels[i]},{row}\n")
