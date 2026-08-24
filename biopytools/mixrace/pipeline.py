"""mixrace 流程步骤(v0.3 GTX 后端)|mixrace pipeline steps (GTX backend).

step01 run_index/run_qc(bwa-mem2 寄主索引由 host_filter 负责;病原索引保留给
fastq2vcf-gtx 复用场景)+ fastp QC;比对与联合 calling 走 gtx_backend;
run_kmer 调 smudgescope;k-mer 谱由 05_kmer/mapped_fastq 输入。
|Steps kept here: genome index + fastp QC + depth cache + smudgescope k-mer.
Mapping/joint calling lives in gtx_backend (fastq2vcf-gtx).
"""
import os
from pathlib import Path
from typing import Optional

from .utils import get_conda_env


def _done(ckpt, step: str, must_exist) -> bool:
    """断点有效 = .done 存在 且 关键输出存在(自愈)。|valid = .done AND key output exists."""
    if not ckpt.exists(step):
        return False
    if must_exist is None:
        return True
    return Path(str(must_exist)).exists()


def run_index(config, runner, ckpt) -> Path:
    """step01a: bwa-mem2 索引参考基因组(全局,一次)|bwa-mem2 index (global, once)."""
    idx_dir = Path(config.output_dir) / "00_pipeline_info" / "index"
    idx_dir.mkdir(parents=True, exist_ok=True)
    fa_copy = idx_dir / os.path.basename(config.genome)
    if config.enable_checkpoint and _done(ckpt, "genome_index", str(fa_copy) + ".0123"):
        runner.logger.info("跳过已完成步骤|Skipping completed step: genome_index")
        return idx_dir
    if not fa_copy.exists():
        import shutil
        shutil.copy(config.genome, str(fa_copy))
    runner.logger.info("开始步骤|Starting step: bwa-mem2 index")
    ok, _, _ = runner.run_conda(config.bwa_mem2_path, ["index", str(fa_copy)],
                                "bwa-mem2索引|bwa-mem2 indexing")
    if config.enable_checkpoint and ok:
        ckpt.create("genome_index")
    elif config.enable_checkpoint:
        runner.logger.error("genome_index 失败,未建断点|genome_index failed, no checkpoint")
    return idx_dir


def run_qc(config, runner, ckpt) -> Optional[Path]:
    """step01b: fastp QC(仅 raw 输入时跑;--clean-fastq-dir 则跳过)|fastp QC (raw only)."""
    if config.clean_fastq_dir:
        runner.logger.info("提供 --clean-fastq-dir,跳过 QC|clean-fastq-dir given, skip QC")
        return Path(config.clean_fastq_dir)
    clean_dir = Path(config.output_dir) / "01_qc"
    if config.enable_checkpoint and _done(ckpt, "qc", clean_dir):
        runner.logger.info("跳过已完成步骤|Skipping completed step: qc")
        return clean_dir
    clean_dir.mkdir(parents=True, exist_ok=True)
    runner.logger.info("开始步骤|Starting step: fastp QC")
    ok, _, _ = runner.run(
        f"biopytools fastp -i {config.fastq_dir} -o {clean_dir} -t {config.threads}",
        "fastp质控|fastp QC")
    if config.enable_checkpoint and ok:
        ckpt.create("qc")
    elif config.enable_checkpoint:
        runner.logger.error("qc 失败,未建断点|qc failed, no checkpoint")
    return clean_dir


def _parse_sn(stats_text: str, key: str) -> Optional[int]:
    """从 samtools stats 取 SN 整数值(字段名精确匹配)|parse SN integer (exact field match).

    子串匹配会先命中 'bases mapped (cigar)' 行,必须整字段相等|substring matching
    hits the 'bases mapped (cigar)' line first; require exact field equality.
    """
    for line in stats_text.splitlines():
        if not line.startswith("SN"):
            continue
        f = line.split("\t")
        if len(f) >= 3 and f[1].rstrip(":") == key:   # SN 字段名带尾冒号|field has trailing ':'
            try:
                return int(f[2])
            except ValueError:
                return None
    return None


def run_depth(runner, config, sample: str, bam: str, genome_size: int):
    """samtools stats → 04_het_eval/alignment_qc/{sample}.stats.txt;平均深度 = bases_mapped / genome_size。
    |samtools stats -> stats.txt; mean depth = bases_mapped / genome_size."""
    qc_dir = Path(config.output_dir) / "04_het_eval" / "alignment_qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    ok, stats_txt, _ = runner.run_conda(config.samtools_path, ["stats", "-@", str(config.threads), str(bam)],
                                        f"samtools stats {sample}")
    if not (ok and stats_txt):
        return None
    (qc_dir / f"{sample}.stats.txt").write_text(stats_txt)
    bases = _parse_sn(stats_txt, "bases mapped")
    if bases and genome_size:
        return bases / genome_size
    return parse_mean_depth(stats_txt)   # 回退:某些 stats 有 average depth 行|fallback


def read_cached_depth(stats_file, genome_size: int):
    """从已缓存的 stats.txt 算平均深度(bases_mapped/genome_size),不重跑 samtools stats。
    |compute mean depth from cached stats.txt (bases_mapped/genome_size), no re-run."""
    try:
        txt = Path(stats_file).read_text(encoding="utf-8")
    except OSError:
        return None
    bases = _parse_sn(txt, "bases mapped")
    if bases and genome_size:
        return bases / genome_size
    return parse_mean_depth(txt)


def run_kmer(config, runner, ckpt, clean_dir: str) -> Path:
    """step04b: k-mer 谱(smudgescope,读 mapped fastq)|k-mer spectrum via smudgescope."""
    kmer_root = Path(config.output_dir) / "05_kmer"
    if config.enable_checkpoint and _done(ckpt, "kmer", kmer_root):
        runner.logger.info("跳过已完成步骤|Skipping completed step: kmer")
        return kmer_root
    runner.logger.info("开始步骤|Starting step: smudgescope k-mer谱")
    ok, _, _ = runner.run(
        f"biopytools smudgescope -i {clean_dir} -o {kmer_root} "
        f"-l {config.read_length} -k {config.kmer_size} -t {config.threads} "
        f"--read1-suffix *_1.mapped.fq.gz",
        "k-mer谱分析|k-mer spectrum (smudgescope)")
    if config.enable_checkpoint and ok:
        ckpt.create("kmer")
    elif config.enable_checkpoint:
        runner.logger.error("kmer 失败,未建断点|kmer failed, no checkpoint")
    return kmer_root


def parse_mean_depth(stats_text: str) -> Optional[float]:
    """从 samtools stats 解析 average depth(回退用)|parse average depth (fallback)."""
    if not stats_text:
        return None
    for line in stats_text.splitlines():
        if line.startswith("SN") and "average depth" in line:
            try:
                return float(line.split("\t")[-1])
            except (ValueError, IndexError):
                return None
    return None


def parse_genomescope_model(text: str) -> dict:
    """解析 genomescope2 model.txt(kcov + 杂合度 r)|parse model.txt."""
    out = {}
    if not text:
        return out
    for line in text.splitlines():
        f = line.replace(",", " ").split()
        if len(f) >= 2:
            k = f[0].lower().rstrip(":")
            try:
                v = float(f[1])
            except ValueError:
                continue
            if k in ("kcov", "kmercov"):
                out["kcov"] = v
            elif k == "r1" or (k == "r" and "heterozygosity" not in out):
                # genomescope2 model.txt: r1 为杂合度(p=2);优先 r1,勿被 r2 覆盖
                out["heterozygosity"] = v
    return out
