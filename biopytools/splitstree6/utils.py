"""splitstree6 工具函数|SplitsTree6 utilities

日志管理器、Xvfb 虚拟显示管理、VCF→距离矩阵转换(p-distance,numpy 向量化)。
|Logger manager, Xvfb virtual display management, VCF→distance conversion
(numpy-vectorized p-distance).
"""

import gzip
import logging
import os
import subprocess
import sys
from typing import Optional, Tuple

import numpy as np


class ModuleLogger:
    """模块日志管理器(三 handler:stdout/stderr/file)|Module logger (3 handlers)"""

    def __init__(self, log_file: Optional[str] = None, log_level: str = "INFO"):
        self.log_file = log_file
        self.logger = logging.getLogger("splitstree6")
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


def read_vcf_gt_matrix(vcf_path: str, logger=None):
    """从 VCF 读取样本基因型并计算 p-distance 矩阵(numpy 向量化)
    |Read VCF genotypes and compute p-distance matrix (numpy vectorized)

    仅处理双等位位点;缺失基因型在配对比较中跳过。
    向量化:每位点按联合等位 (allele1, allele2) 编码为单个整数 3*a1+a2
    (取值 0..8),p-distance = 联合等位不同的位点数 / 共享位点数。
    |Each site encoded as combined genotype 3*a1+a2 (0..8); a site differs
    when the combined genotypes differ (equivalent to GT string comparison).

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

    # 共享位点数|shared count per pair
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
    """写 SplitsTree6 兼容的距离 CSV(不带首行样本数)|Write SplitsTree6-compatible distance CSV

    关键:SplitsTree6 的 CSVReader.acceptsFirstLine 要求**首行**含 ≥4 个逗号,
    纯数字的首行(样本数)会被判为未知格式 → 必须省略。
    |SplitsTree6's CSVReader.acceptsFirstLine requires commas on the FIRST line;
    a leading sample-count line breaks auto-detection — omit it.
    """
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    n = len(labels)
    with open(out_path, "w") as fh:
        for i in range(n):
            row = ",".join(f"{v:.6f}" for v in matrix[i])
            fh.write(f"{labels[i]},{row}\n")


class XvfbDisplay:
    """Xvfb 虚拟显示管理器|Xvfb virtual display manager

    SplitsTree6 主类依赖 JavaFX GUI 栈,须有 X DISPLAY 才能启动;
    用 Xvfb 提供虚拟显示实现 headless 运行。DISPLAY 已存在时直接复用。
    |SplitsTree6 requires JavaFX (a display); run under Xvfb headlessly,
    reusing an existing display when present.
    """

    def __init__(self, xvfb_path: str, display: str = ":99",
                 screen: str = "1024x768x24", logger=None,
                 lib_hints: Optional[list] = None):
        self.xvfb_path = xvfb_path
        self.display = display
        self.screen = screen
        self.logger = logger
        self.lib_hints = [os.path.expanduser(p) for p in (lib_hints or [])]
        self.proc = None

    def _env(self) -> dict:
        """带 LD_LIBRARY_PATH 提示的环境|Environment with lib hints prepended"""
        env = dict(os.environ)
        hints = os.pathsep.join(
            p for p in [*self.lib_hints, env.get("LD_LIBRARY_PATH", "")] if p)
        if hints:
            env["LD_LIBRARY_PATH"] = hints
        return env

    def start(self) -> bool:
        """启动或复用 Xvfb|Start Xvfb or reuse an existing one"""
        cmd = [self.xvfb_path, self.display, "-screen", "0", self.screen]
        if self.logger:
            self.logger.info(f"命令|Command: {' '.join(cmd)}")
        try:
            self.proc = subprocess.Popen(
                cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                env=self._env())
        except FileNotFoundError:
            if self.logger:
                self.logger.error(f"Xvfb 未找到|Xvfb not found: {self.xvfb_path}")
            return False
        time.sleep(2.0)   # 等 Xvfb 完成初始化|let Xvfb finish initializing
        return True       # 显示不可用时由后续真实运行报错|failures surface downstream

    @property
    def env(self) -> dict:
        """带 DISPLAY 的环境变量|Environment dict with DISPLAY set"""
        env = dict(os.environ)
        env["DISPLAY"] = self.display
        return env

    def stop(self) -> None:
        """停止本次启动的 Xvfb|Stop the Xvfb we started (if any)"""
        if self.proc and self.proc.poll() is None:
            self.proc.terminate()
            try:
                self.proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.proc.kill()
