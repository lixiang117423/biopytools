"""xpclr 核心计算|xpclr calculator.

VCF header 解析 → 样本校验 → 逐染色体调用 xpclr CLI(断点续传+容错)
→ 合并全基因组表 → Top 候选窗口汇总 → software_versions.yml。
|Parse VCF header → validate samples → per-chrom xpclr CLI (checkpoint +
graceful degradation) → genome-wide merge → top windows → versions yml.
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

from ..common.conda_runner import CommandRunner, build_conda_command
from . import __version__
from .config import XpclrConfig
from .utils import detect_xpclr_version, format_number

# ##contig=<ID=xxx,...> 提取 ID|extract contig ID
_CONTIG_RE = re.compile(r"##contig=.*?<ID=([^,>]+)")


def read_sample_list(path: str) -> List[str]:
    """读样本列表(每行一个ID,跳空行)|Read sample IDs (one per line, blanks skipped)."""
    samples: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if s:
                samples.append(s)
    return samples


def parse_vcf_header(vcf_path: str) -> Tuple[List[str], List[str]]:
    """流式解析 VCF header → (contig 列表, 样本列表)|Stream header → (contigs, samples).

    遇第一条非 # 行即停(不读数据区);contig 按 header 出现顺序;
    .gz 用 gzip 读取,普通文本 VCF 直接打开(防御性容错,不依赖调用方保证)。
    |Stops at first record line; contigs in header order; gzip for .gz,
    plain open for uncompressed VCF (defensive).
    """
    contigs: List[str] = []
    samples: List[str] = []
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##contig="):
                m = _CONTIG_RE.search(line)
                if m:
                    contigs.append(m.group(1))
            elif line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
    return contigs, samples


def validate_samples(samples_a: List[str], samples_b: List[str],
                     vcf_samples: List[str]) -> None:
    """A/B 重叠与幽灵样本,一次性收集全部错误|Collect overlap/ghost errors at once."""
    errors: List[str] = []
    overlap = sorted(set(samples_a) & set(samples_b))
    if overlap:
        errors.append(
            f"A/B 样本重叠(结果无意义)|A/B sample overlap: {','.join(overlap)}")
    vcf_set = set(vcf_samples)
    for name, group in (("A", samples_a), ("B", samples_b)):
        missing = sorted({s for s in group if s not in vcf_set})
        if missing:
            errors.append(
                f"群体{name}样本不在VCF中|Pop {name} samples missing from VCF: "
                f"{','.join(missing)}")
    if errors:
        raise ValueError("\n".join(errors))


class XpclrCalculator:
    """xpclr 编排器|xpclr orchestrator."""

    def __init__(self, config: XpclrConfig, logger):
        self.config = config
        self.logger = logger

    def resolve_chroms(self, contigs: List[str]) -> List[str]:
        """--chroms 过滤(保持 VCF header 顺序);缺省全部|Filter by --chroms in VCF order."""
        if not self.config.chroms:
            return list(contigs)
        contig_set = set(contigs)
        unknown = [c for c in self.config.chroms if c not in contig_set]
        if unknown:
            raise ValueError(
                f"--chroms 含 VCF 不存在的染色体|chroms not in VCF: {','.join(unknown)}")
        order = {c: i for i, c in enumerate(contigs)}
        return sorted(self.config.chroms, key=lambda c: order[c])

    # -- 逐染色体执行|per-chromosome execution --------------------------

    def chrom_output(self, chrom: str) -> Path:
        """单染色体输出路径|Per-chrom output path."""
        return Path(self.config.output_dir) / "01_xpclr" / f"{chrom}.xpclr.tsv"

    def _build_args(self, chrom: str, out_tsv: Path) -> List[str]:
        """拼 xpclr CLI 参数|Build xpclr CLI args."""
        args = [
            "--format", "vcf",
            "--input", self.config.input_vcf,
            "--samplesA", self.config.samples_a,
            "--samplesB", self.config.samples_b,
            "--chr", chrom,
            "--out", str(out_tsv),
            "--size", str(self.config.size),
            "--step", str(self.config.step),
            "--maxsnps", str(self.config.maxsnps),
            "--minsnps", str(self.config.minsnps),
            "--ld", str(self.config.ld),
            "--rrate", str(self.config.rrate),
        ]
        if self.config.phased:
            args.append("--phased")
        return args

    def run_chrom(self, chrom: str, runner: CommandRunner) -> bool:
        """单染色体:断点续传 + 失败 WARNING 继续|Checkpoint + warn-and-continue."""
        out_tsv = self.chrom_output(chrom)
        if out_tsv.exists() and out_tsv.stat().st_size > 0:
            self.logger.info(f"跳过已完成染色体|Skipping completed chrom: {chrom}")
            return True
        cmd = build_conda_command(self.config.xpclr_path,
                                  self._build_args(chrom, out_tsv))
        ok, _, stderr = runner.run(cmd, f"XP-CLR 染色体|XP-CLR chrom {chrom}")
        if not ok or not out_tsv.exists() or out_tsv.stat().st_size == 0:
            # 失败清理残留半成品,防止重跑被"非空即完成"误跳过丢数据
            # |drop partial output on failure so the checkpoint check won't skip it
            try:
                out_tsv.unlink()
            except OSError:
                pass
            self.logger.warning(
                f"染色体 {chrom} 失败,跳过继续|Chrom {chrom} failed, skipped. "
                f"stderr: {(stderr or '')[:300]}")
            return False
        return True

    # -- 合并与汇总|merge & summary --------------------------------------

    def merge_tables(self, chroms: List[str]) -> Optional[Path]:
        """合并成功染色体的 TSV|Merge per-chrom TSVs into genome-wide table."""
        out_dir = Path(self.config.output_dir) / "02_merged"
        out_dir.mkdir(parents=True, exist_ok=True)
        frames = []
        for chrom in chroms:
            try:
                frames.append(pd.read_csv(self.chrom_output(chrom), sep="\t"))
            except pd.errors.EmptyDataError:
                self.logger.warning(
                    f"染色体 {chrom} 输出为空,合并时跳过|Chrom {chrom} output empty, "
                    f"skipped in merge")
        if not frames:
            return None
        merged = pd.concat(frames, ignore_index=True)
        out = out_dir / f"{self.config.label}.xpclr.genome.tsv"
        merged.to_csv(out, sep="\t", index=False)
        self.logger.info(
            f"全基因组合并表|Genome-wide merged table: {out} "
            f"(窗口数|windows: {format_number(len(merged))})")
        return out

    def top_windows(self, merged_path: Path) -> Optional[Path]:
        """按 xpclr_norm 降序取 Top N|Top N windows by xpclr_norm descending."""
        df = pd.read_csv(merged_path, sep="\t")
        if "xpclr_norm" not in df.columns or df["xpclr_norm"].isna().all():
            self.logger.warning(
                "无有效 xpclr_norm,跳过 Top 汇总|No valid xpclr_norm, top summary skipped")
            return None
        top = (df.dropna(subset=["xpclr_norm"])
                 .sort_values("xpclr_norm", ascending=False)
                 .head(self.config.top_n))
        out_dir = Path(self.config.output_dir) / "03_top"
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / f"{self.config.label}.xpclr.top{self.config.top_n}.tsv"
        top.to_csv(out, sep="\t", index=False)
        self.logger.info(
            f"Top 候选窗口表|Top candidate windows: {out} (n={len(top)})")
        return out

    def write_versions(self) -> None:
        """写 software_versions.yml|Write software_versions.yml."""
        import yaml
        out_dir = Path(self.config.output_dir) / "00_pipeline_info"
        out_dir.mkdir(parents=True, exist_ok=True)
        info = {
            "biopytools_xpclr": {"version": __version__},
            "xpclr": {"version": detect_xpclr_version(self.config.xpclr_path),
                      "path": self.config.xpclr_path},
            "parameters": {
                "input_vcf": self.config.input_vcf,
                "label": self.config.label,
                "size": self.config.size,
                "step": self.config.step,
                "maxsnps": self.config.maxsnps,
                "minsnps": self.config.minsnps,
                "ld": self.config.ld,
                "phased": self.config.phased,
                "rrate": self.config.rrate,
                "top_n": self.config.top_n,
            },
        }
        (out_dir / "software_versions.yml").write_text(
            yaml.safe_dump(info, allow_unicode=True, sort_keys=False), encoding="utf-8")

    # -- 总编排|orchestration --------------------------------------------

    def run(self) -> bool:
        """端到端编排:校验→逐染色体→合并→Top→版本|End-to-end orchestration."""
        contigs, vcf_samples = parse_vcf_header(self.config.input_vcf)
        if not contigs:
            raise ValueError("VCF header 无 ##contig 行|No ##contig lines in VCF header")
        samples_a = read_sample_list(self.config.samples_a)
        samples_b = read_sample_list(self.config.samples_b)
        if not samples_a or not samples_b:
            raise ValueError(
                f"样本列表为空|Sample list empty: "
                f"A={len(samples_a)} B={len(samples_b)}")
        validate_samples(samples_a, samples_b, vcf_samples)
        chroms = self.resolve_chroms(contigs)
        self.logger.info(
            f"染色体数|Chromosomes: {len(chroms)}; "
            f"样本|Samples: A={len(samples_a)} B={len(samples_b)}")
        Path(self.config.output_dir, "01_xpclr").mkdir(parents=True, exist_ok=True)
        runner = CommandRunner(self.logger)
        ok_chroms = [c for c in chroms if self.run_chrom(c, runner)]
        failed = len(chroms) - len(ok_chroms)
        self.logger.info(
            f"染色体完成|Chroms done: {len(ok_chroms)} 成功|ok, {failed} 失败|failed")
        if not ok_chroms:
            self.logger.error("全部染色体失败,中止|All chromosomes failed, aborting")
            return False
        merged = self.merge_tables(ok_chroms)
        if merged is None:
            self.logger.error("无可合并结果,中止|Nothing to merge, aborting")
            return False
        top = self.top_windows(merged)
        if top:
            self.logger.info(f"Top 候选窗口表|Top windows: {top}")
        self.write_versions()
        if failed:
            self.logger.warning(
                f"{failed} 条染色体失败,结果为部分基因组|{failed} chrom(s) failed; "
                f"results cover a partial genome")
        return True
