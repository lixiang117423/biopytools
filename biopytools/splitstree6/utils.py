"""splitstree6 工具函数|SplitsTree6 utilities

日志管理器、Xvfb 虚拟显示管理、VCF→距离矩阵转换(p-distance)。
|Logger manager, Xvfb virtual display management, VCF→distance conversion.
"""

import gzip
import logging
import os
import subprocess
import sys
import time
from typing import Optional, Tuple


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


def read_vcf_gt_matrix(vcf_path: str):
    """从 VCF 读取样本基因型并计算 p-distance 矩阵|Read VCF genotypes and compute p-distances

    仅处理双等位位点;缺失基因型在配对比较中跳过。
    |Biallelic sites only; missing genotypes skipped per-pair.

    Returns:
        (labels, matrix): 样本名列表与 n×n 距离矩阵|labels and n×n distance matrix
    """
    labels = []
    genotypes = []
    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                labels = parts[9:]
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            gts = []
            for field in parts[9:]:
                gt = field.split(":")[0]
                alleles = gt.replace("|", "/").split("/")
                if len(alleles) < 2 or "." in alleles:
                    gts.append(None)
                else:
                    try:
                        gts.append((int(alleles[0]), int(alleles[1])))
                    except ValueError:
                        gts.append(None)
            genotypes.append(gts)

    n = len(labels)
    matrix = [[0.0] * n for _ in range(n)]
    if n == 0:
        raise ValueError("VCF 中未找到样本|No samples found in VCF header")

    # 配对计算 p-distance|pairwise p-distance computation
    for i in range(n):
        matrix[i][i] = 0.0
        for j in range(i + 1, n):
            diff = same = 0
            for gts in genotypes:
                a, b = gts[i], gts[j]
                if a is None or b is None:
                    continue
                same += 1
                if a != b:
                    diff += 1
            dist = (diff / same) if same else float("nan")
            matrix[i][j] = dist
            matrix[j][i] = dist

    # NaN(无共享位点)置 0|NaN (no shared sites) → 0
    for i in range(n):
        for j in range(n):
            if matrix[i][j] != matrix[i][j]:   # NaN check
                matrix[i][j] = 0.0
    return labels, matrix


def write_distance_csv(labels, matrix, out_path: str) -> None:
    """写 CSV 变体(a):首行样本数,随后 label + n 个值——SplitsTree6 CSVReader 兼容
    |Write CSV variant (a): first line sample count, then label + values per row
    """
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    n = len(labels)
    with open(out_path, "w") as fh:
        fh.write(f"{n}\n")
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
