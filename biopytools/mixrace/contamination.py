"""mixrace 污染评估(step 6)|mixrace contamination assessment (kraken2 + bracken).

kraken2 逐样本分类 clean reads:仅留 --report,--output 指向 /dev/null——
kraken2 默认把逐 read 标签打到 stdout,而 CommandRunner 捕获 stdout,
深度样本会有 GB 级文本灌爆管道。bracken 按层级重估丰度(-r 吸附到库内
可用 kmer_distrib 读长)。汇总为 contamination_summary.tsv / detail.tsv,
独立于 verdict_summary,不触碰既有输出。样本间串行:DB 每进程整载,
并行 = 内存翻 N 倍。
|kraken2 per-sample classification of clean reads: --report only, --output
/dev/null — kraken2 prints per-read labels to stdout by default, which would
flood CommandRunner's captured pipe on deep samples. bracken re-estimates
abundance per rank (-r snapped to an available kmer_distrib length).
Standalone summary tables; verdict_summary untouched. Samples sequential:
the DB is loaded per process, parallelism multiplies RAM.
"""
import os
from pathlib import Path
from typing import Dict, List, Optional

# 数据库自带 kmer_distrib 读长,bracken -r 必须命中其一
# |read lengths shipped with the DB; bracken -r must match one of them
KMER_DISTRIB_LENGTHS = (50, 75, 100, 150, 200, 250, 300)


def snap_read_length(read_length: int, available=KMER_DISTRIB_LENGTHS) -> int:
    """吸附到最近可用 kmer_distrib 读长|snap to nearest available read length."""
    return min(available, key=lambda n: (abs(n - read_length), n))


def parse_kraken2_report(report_file) -> Dict[str, int]:
    """解析 kraken2 report 的 U(unclassified)/R(root) 行|parse U/root rows.

    返回 {"total_reads","classified_reads","unclassified_reads"};
    root 行缺失(空 report)或文件不存在 → {}。
    |returns those counts; {} on missing root row or missing file.
    """
    try:
        text = Path(report_file).read_text(encoding="utf-8")
    except OSError:
        return {}
    classified = None
    unclassified = 0
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        if f[3] == "R" and f[4] == "1":
            try:
                classified = int(f[1])
            except ValueError:
                continue
        elif f[3] == "U":
            try:
                unclassified = int(f[1])
            except ValueError:
                continue
    if classified is None:
        return {}
    return {"total_reads": classified + unclassified,
            "classified_reads": classified,
            "unclassified_reads": unclassified}


def parse_bracken_report(report_file) -> List[Dict]:
    """解析 bracken 丰度表(7 列)|parse bracken abundance table (7 columns)."""
    try:
        lines = Path(report_file).read_text(encoding="utf-8").splitlines()
    except OSError:
        return []
    rows: List[Dict] = []
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) < 7:
            continue
        try:
            rows.append({
                "name": f[0], "taxonomy_id": f[1], "taxonomy_lvl": f[2],
                "kraken_assigned_reads": int(f[3]), "added_reads": int(f[4]),
                "new_est_reads": int(f[5]),
                "fraction_total_reads": float(f[6]),
            })
        except ValueError:
            continue
    return rows


def _pct(x: float) -> str:
    return f"{x * 100:.2f}"


def render_summary_tsv(per_sample: Dict[str, Dict]) -> str:
    """每样本一行汇总表|one row per sample summary table.

    物种占比自算 new_est_reads/k2总reads(含unclassified)。bracken 原生
    fraction_total_reads 分母只有 classified reads(est_abundance.py 的
    sum_all_reads),与本表其余列口径不一致,不采信。
    |species pct recomputed as new_est_reads / k2 total (incl. unclassified);
    bracken's native fraction uses a classified-only denominator.
    """
    header = ["sample", "total_reads", "classified_pct", "unclassified_pct",
              "top_species", "top_species_pct", "other_classified_pct",
              "n_species_ge_1pct"]
    lines = ["\t".join(header)]
    for sample in sorted(per_sample):
        k2 = per_sample[sample].get("k2") or {}
        bracken = per_sample[sample].get("bracken") or []
        if not k2:
            lines.append("\t".join([sample] + ["NA"] * 7))
            continue
        total = k2["total_reads"]
        classified_pct = k2["classified_reads"] / total if total else 0.0
        unclassified_pct = k2["unclassified_reads"] / total if total else 0.0
        if bracken and total:
            top = max(bracken, key=lambda r: r["new_est_reads"])
            top_name = top["name"]
            top_pct = top["new_est_reads"] / total
        else:
            top_name, top_pct = "NA", 0.0
        other_pct = max(0.0, classified_pct - top_pct)
        n_ge1 = sum(1 for r in bracken
                    if total and r["new_est_reads"] / total >= 0.01)
        lines.append("\t".join([
            sample, str(total), _pct(classified_pct), _pct(unclassified_pct),
            top_name, _pct(top_pct), _pct(other_pct), str(n_ge1)]))
    return "\n".join(lines) + "\n"


def render_detail_tsv(per_sample: Dict[str, Dict], min_fraction: float = 0.001) -> str:
    """长表:样本×物种(占比≥min_fraction)|long table: sample x species.

    fraction_total_reads 列 = new_est_reads/k2总reads(含unclassified),
    口径与汇总表一致;k2 缺失的样本无法归一,跳过。
    |fraction column = new_est_reads / k2 total (incl. unclassified),
    consistent with the summary; samples without k2 counts are skipped.
    """
    header = ["sample", "name", "taxonomy_id", "taxonomy_lvl",
              "new_est_reads", "fraction_total_reads"]
    lines = ["\t".join(header)]
    for sample in sorted(per_sample):
        k2 = per_sample[sample].get("k2") or {}
        total = k2.get("total_reads") or 0
        if not total:
            continue
        for r in per_sample[sample].get("bracken") or []:
            frac = r["new_est_reads"] / total
            if frac < min_fraction:
                continue
            lines.append("\t".join([
                sample, r["name"], r["taxonomy_id"], r["taxonomy_lvl"],
                str(r["new_est_reads"]), f"{frac:.6f}"]))
    return "\n".join(lines) + "\n"


def _warn_ram_needs(runner, config) -> None:
    """提示 kraken2 整库载入内存的实际需求|log actual full-load RAM need."""
    if config.kraken_memory_mapping:
        runner.logger.info("kraken2 --memory-mapping 已启用(省内存,较慢)|"
                           "kraken2 memory-mapping enabled (less RAM, slower)")
        return
    size_note = ""
    try:
        gb = os.path.getsize(Path(config.kraken2_db) / "hash.k2d") / 1e9
        size_note = f"(hash.k2d {gb:.0f}GB)"
    except OSError:
        pass
    runner.logger.warning(
        f"kraken2 默认将数据库整载入内存{size_note},确认节点内存充足,"
        f"不足时用 --kraken-memory-mapping|kraken2 loads the full DB into RAM "
        f"{size_note}; use --kraken-memory-mapping on smaller nodes")


def run_kraken2_sample(config, runner, ckpt, sample: str,
                       r1: str, r2: str) -> Optional[Path]:
    """kraken2 单样本分类|kraken2 classification for one sample."""
    out_dir = Path(config.output_dir) / "08_contamination" / "kraken2"
    out_dir.mkdir(parents=True, exist_ok=True)
    report = out_dir / f"{sample}_k2_report.txt"
    step = f"kraken2_{sample}"
    if config.enable_checkpoint and ckpt.exists(step) and report.exists():
        runner.logger.info(f"跳过已完成步骤|Skipping completed step: {step}")
        return report
    args = ["--db", config.kraken2_db, "--paired",
            "--threads", str(config.threads),
            "--report", str(report),
            "--output", "/dev/null"]
    if config.kraken_memory_mapping:
        args.append("--memory-mapping")
    args += [r1, r2]
    runner.logger.info(f"开始步骤|Starting step: kraken2 {sample}")
    ok, _, _ = runner.run_conda(
        config.kraken2_path, args, f"kraken2分类|kraken2 classification {sample}")
    if ok and report.exists():
        if config.enable_checkpoint:
            ckpt.create(step)
    else:
        runner.logger.error(f"kraken2 失败,样本进入汇总为NA|kraken2 failed: {sample}")
    return report if (ok and report.exists()) else None


def run_bracken_sample(config, runner, ckpt, sample: str,
                       k2_report: Path) -> Optional[Path]:
    """bracken 单样本丰度重估|bracken abundance re-estimation for one sample."""
    out_dir = Path(config.output_dir) / "08_contamination" / "bracken"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{sample}_bracken_{config.bracken_level}.txt"
    step = f"bracken_{sample}"
    if config.enable_checkpoint and ckpt.exists(step) and out.exists():
        runner.logger.info(f"跳过已完成步骤|Skipping completed step: {step}")
        return out
    read_len = snap_read_length(config.read_length)
    if read_len != config.read_length:
        runner.logger.warning(
            f"read_length={config.read_length} 无对应 kmer_distrib,"
            f"bracken -r 吸附到|snapped to: {read_len}")
    args = ["-d", config.kraken2_db, "-i", str(k2_report), "-o", str(out),
            "-r", str(read_len), "-t", str(config.threads),
            "-l", config.bracken_level]
    runner.logger.info(f"开始步骤|Starting step: bracken {sample}")
    ok, _, _ = runner.run_conda(
        config.bracken_path, args, f"bracken丰度|bracken abundance {sample}")
    if ok and out.exists():
        if config.enable_checkpoint:
            ckpt.create(step)
    else:
        runner.logger.error(f"bracken 失败,样本进入汇总为NA|bracken failed: {sample}")
    return out if (ok and out.exists()) else None


def run_contamination(config, runner, ckpt, samples) -> None:
    """step6 编排:逐样本 kraken2→bracken→汇总|orchestrate per sample + summary."""
    if not samples:
        runner.logger.warning("污染评估无样本可跑|no samples for contamination assessment")
        return
    contam_dir = Path(config.output_dir) / "08_contamination"
    contam_dir.mkdir(parents=True, exist_ok=True)
    _warn_ram_needs(runner, config)
    per_sample: Dict[str, Dict] = {}
    for s in samples:
        sample = s["sample"]
        k2_report = run_kraken2_sample(config, runner, ckpt, sample,
                                       s["r1"], s["r2"])
        if k2_report is None:
            per_sample[sample] = {}
            continue
        bracken_out = run_bracken_sample(config, runner, ckpt, sample, k2_report)
        per_sample[sample] = {
            "k2": parse_kraken2_report(k2_report),
            "bracken": parse_bracken_report(bracken_out) if bracken_out else [],
        }
    (contam_dir / "contamination_summary.tsv").write_text(
        render_summary_tsv(per_sample), encoding="utf-8")
    (contam_dir / "contamination_detail.tsv").write_text(
        render_detail_tsv(per_sample), encoding="utf-8")
    runner.logger.info(f"污染评估汇总已写|contamination summary written: "
                       f"{contam_dir / 'contamination_summary.tsv'}")
