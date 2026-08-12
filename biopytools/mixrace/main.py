"""mixrace 主入口与编排|mixrace main entry & orchestration.

解析参数 → 配置 → 按 --step 调度 8 步:01 索引+QC / 02 比对 / 03 变异 / 04 过滤 /
05 VAF / 06 多等位+Fws / 07 k-mer(smudgescope) / 08 判读+报告。逐样本循环 step02-06,
队列 Fws,可选 --pure-samples 校准阈值,输出跨样品对比表。
|argparse -> config -> --step dispatch over 8 steps. Per-sample loop for 02-06,
cohort Fws, optional --pure-samples calibration, cross-sample summary table.
"""
import argparse
import os
import sys
from pathlib import Path

from .config import MixraceConfig
from .utils import ModuleLogger, CommandRunner, CheckpointManager, write_software_versions
from .samples import discover_samples
from .pipeline import (run_index, run_qc, run_align, run_call, run_filter, run_kmer,
                       parse_mean_depth, parse_genomescope_model)
from .vaf_analysis import parse_vcf_ad, compute_vaf, vaf_metrics, multiallelic_stats, compute_fws
from .verdict import judge, calibrate_thresholds
from .reporter import (generate_vaf_histogram_r, build_sample_report, build_summary_table)


def _argv_to_config() -> MixraceConfig:
    """解析命令行参数 → MixraceConfig|parse argv into MixraceConfig."""
    p = argparse.ArgumentParser(prog="mixrace",
                                description="WGS混合小种检测|WGS mixed-race detection")
    p.add_argument("-i", "--input", dest="fastq_dir", required=True)
    p.add_argument("-g", "--genome", required=True)
    p.add_argument("-o", "--output-dir", dest="output_dir", default="mixrace_out")
    p.add_argument("--repeat-bed", dest="repeat_bed", default=None)
    p.add_argument("-t", "--threads", type=int, default=12)
    p.add_argument("-k", "--kmer-size", type=int, default=21)
    p.add_argument("-l", "--read-length", type=int, default=150)
    p.add_argument("--step", type=int, default=None)
    p.add_argument("--no-checkpoint", dest="enable_checkpoint", action="store_false")
    p.add_argument("--dry-run", dest="dry_run", action="store_true")
    p.add_argument("--min-qual", type=int, default=30)
    p.add_argument("--min-dp", type=int, default=15)
    p.add_argument("--min-alt-reads", type=int, default=3)
    p.add_argument("--vaf-mid-low", type=float, default=0.05)
    p.add_argument("--vaf-mid-high", type=float, default=0.15)
    p.add_argument("--multiallelic-low", type=float, default=0.005)
    p.add_argument("--multiallelic-high", type=float, default=0.03)
    p.add_argument("--fws-cutoff", type=float, default=0.95)
    p.add_argument("--min-depth", type=float, default=15.0)
    p.add_argument("--pure-samples", dest="pure_samples", default=None)
    a = p.parse_args()
    pure = a.pure_samples.split(",") if a.pure_samples else None
    return MixraceConfig(
        fastq_dir=a.fastq_dir, genome=a.genome, output_dir=a.output_dir,
        repeat_bed=a.repeat_bed, threads=a.threads, kmer_size=a.kmer_size,
        read_length=a.read_length, step=a.step, enable_checkpoint=a.enable_checkpoint,
        dry_run=a.dry_run, min_qual=a.min_qual, min_dp=a.min_dp,
        min_alt_reads=a.min_alt_reads, vaf_mid_low=a.vaf_mid_low,
        vaf_mid_high=a.vaf_mid_high, multiallelic_low=a.multiallelic_low,
        multiallelic_high=a.multiallelic_high, fws_cutoff=a.fws_cutoff,
        min_depth=a.min_depth, pure_samples=pure)


def _sample_vaf_metrics(config, runner, sample: str, filt_vcf: str):
    """从 filtered VCF 提取 AD → VAF/Hw/多等位;写 vaf.tsv;返回 (metrics_partial, per_site_vafs)。
    |extract AD from filtered VCF -> VAF/Hw/multiallelic; write vaf.tsv; return metrics + per-site."""
    ok, qtxt, _ = runner.run_conda(
        config.bcftools_path,
        ["query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n", str(filt_vcf)],
        f"query AD|query AD {sample}")
    recs = parse_vcf_ad(qtxt)
    vaf_recs = compute_vaf(recs)
    vaf_dir = Path(config.output_dir) / sample / "05_vaf"
    vaf_dir.mkdir(parents=True, exist_ok=True)
    vaf_tsv = vaf_dir / f"{sample}.vaf.tsv"
    with open(vaf_tsv, "w") as fh:
        fh.write("chrom\tpos\tref\talts\tdp\tad_ref\tad_alts\tvafs\tis_multi\n")
        for r in vaf_recs:
            fh.write("\t".join([
                r["chrom"], str(r["pos"]), r["ref"], ",".join(r["alts"]),
                str(r["dp"]), str(r["ad_ref"]), ",".join(str(x) for x in r["ad_alts"]),
                ",".join(f"{v:.4f}" for v in r["vafs"]),
                "1" if r["is_multiallelic"] else "0",
            ]) + "\n")
    all_vafs = [v for r in vaf_recs for v in r["vafs"]]
    vm = vaf_metrics(all_vafs)
    ma = multiallelic_stats(recs)
    per_site = {(r["chrom"], r["pos"]): (r["vafs"][0] if r["vafs"] else 0.0) for r in vaf_recs}
    return ({"vaf_mid_ratio": vm["intermediate_ratio"], "hw": vm["hw"],
             "multiallelic_ratio": ma["ratio"]}, per_site)


def _read_mean_depth(config, sample: str):
    """从 02_alignment/stats.txt 读平均深度|read mean depth from stats.txt."""
    stats = Path(config.output_dir) / sample / "02_alignment" / "stats.txt"
    if stats.exists():
        return parse_mean_depth(stats.read_text())
    return None


def _read_heterozygosity(config, sample: str):
    """从 07_kmer/<sample>/02_genomescope/model.txt 读杂合度|read het from model.txt."""
    model = Path(config.output_dir) / "07_kmer" / sample / "02_genomescope" / "model.txt"
    if model.exists():
        return parse_genomescope_model(model.read_text()).get("heterozygosity")
    return None


def run_pipeline(config, runner, ckpt, logger):
    """8 步编排|orchestrate 8 steps by --step."""
    step = config.step
    thr = {"vaf_mid_low": config.vaf_mid_low, "vaf_mid_high": config.vaf_mid_high,
           "multiallelic_low": config.multiallelic_low, "multiallelic_high": config.multiallelic_high,
           "fws_cutoff": config.fws_cutoff, "min_depth": config.min_depth}

    # 01 索引 + QC|index + QC
    if step in (None, 1):
        idx_dir = run_index(config, runner, ckpt)
        clean_dir = run_qc(config, runner, ckpt)
    else:
        idx_dir = Path(config.output_dir) / "00_pipeline_info" / "index"
        clean_dir = Path(config.output_dir) / "01_qc"

    # 样本列表(优先 clean 目录)|sample list (prefer clean dir)
    samples = []
    if step in (None, 2, 3, 4, 5, 6, 8):
        src = str(clean_dir) if Path(str(clean_dir)).is_dir() else config.fastq_dir
        samples = discover_samples(src)

    fa = str(Path(str(idx_dir)) / os.path.basename(config.genome))
    per_site_all = {}
    partial = {}
    # 02-06 逐样本|per-sample
    for s in samples:
        sample = s["sample"]
        if step in (None, 2):
            md_bam = run_align(config, runner, ckpt, sample, s["r1"], s["r2"], str(idx_dir))
        else:
            md_bam = Path(config.output_dir) / sample / "02_alignment" / f"{sample}.sorted.markdup.bam"
        if step in (None, 3):
            raw = run_call(config, runner, ckpt, sample, str(md_bam), fa)
        else:
            raw = Path(config.output_dir) / sample / "03_variants" / f"{sample}.raw.vcf.gz"
        if step in (None, 4):
            filt = run_filter(config, runner, ckpt, sample, str(raw))
        else:
            filt = Path(config.output_dir) / sample / "04_filtered" / f"{sample}.filtered.vcf.gz"
        if step in (None, 5, 6):
            mp, per_site = _sample_vaf_metrics(config, runner, sample, str(filt))
            partial[sample] = mp
            per_site_all[sample] = per_site
            if step in (None, 5):
                vaf_tsv = Path(config.output_dir) / sample / "05_vaf" / f"{sample}.vaf.tsv"
                png = Path(config.output_dir) / sample / "05_vaf" / f"{sample}.vaf_histogram.png"
                script = generate_vaf_histogram_r(str(vaf_tsv), str(png), config.rscript_path)
                runner.run_conda(config.rscript_path, [script], f"VAF直方图|VAF histogram {sample}")

    # 06 Fws(队列)|cohort Fws
    fws_map = {}
    if step in (None, 5, 6) and len(per_site_all) >= 2:
        fws_map = compute_fws(per_site_all)

    # 07 k-mer 谱|k-mer spectrum (smudgescope)
    if step in (None, 7):
        run_kmer(config, runner, ckpt, str(clean_dir))

    # 08 判读 + 报告|verdict + report
    if step in (None, 8):
        rows = _verdict_and_report(config, runner, logger, samples, partial, per_site_all,
                                   fws_map, thr)
        summ_dir = Path(config.output_dir) / "summary"
        summ_dir.mkdir(parents=True, exist_ok=True)
        tsv, html = build_summary_table(rows)
        (summ_dir / "verdict_summary.tsv").write_text(tsv)
        (summ_dir / "verdict_summary.html").write_text(html)
        logger.info(f"汇总表已写|summary written: {summ_dir / 'verdict_summary.tsv'}")

    # 版本信息|software versions
    write_software_versions(
        config, logger,
        str(Path(config.output_dir) / "00_pipeline_info" / "software_versions.yml"))


def _verdict_and_report(config, runner, logger, samples, partial, per_site_all, fws_map, thr):
    """step08: 逐样品判读 + 报告(支持单独重跑:缺指标则从已有 VCF 重算)|step08 verdict+report.
    Supports standalone rerun: recomputes metrics from existing filtered VCF if missing."""
    # 单独跑 step8 时补算缺失指标|recompute missing metrics for standalone step8
    for s in samples:
        sample = s["sample"]
        if sample not in partial:
            filt = Path(config.output_dir) / sample / "04_filtered" / f"{sample}.filtered.vcf.gz"
            if filt.exists():
                mp, per_site = _sample_vaf_metrics(config, runner, sample, str(filt))
                partial[sample] = mp
                per_site_all[sample] = per_site
    if not fws_map and len(per_site_all) >= 2:
        fws_map = compute_fws(per_site_all)

    # 校准(已知纯样品 mean+2SD)|calibration via known-pure samples
    if config.pure_samples:
        pure_rows = [partial[p] for p in config.pure_samples if p in partial]
        if len(pure_rows) >= 2:
            cal = calibrate_thresholds(pure_rows)
            if cal.get("vaf_mid_low") is not None:
                thr["vaf_mid_low"] = cal["vaf_mid_low"]
            if cal.get("multiallelic_low") is not None:
                thr["multiallelic_low"] = cal["multiallelic_low"]
            logger.info(f"阈值已用 {len(pure_rows)} 个纯样品校准|thresholds calibrated with {len(pure_rows)} pure samples")

    rows = []
    for s in samples:
        sample = s["sample"]
        mp = partial.get(sample, {})
        depth = _read_mean_depth(config, sample)
        het = _read_heterozygosity(config, sample)
        fw = fws_map.get(sample, {}).get("fws")
        hw = fws_map.get(sample, {}).get("hw", mp.get("hw", 0.0))
        metrics = {
            "vaf_mid_ratio": mp.get("vaf_mid_ratio", 0.0),
            "multiallelic_ratio": mp.get("multiallelic_ratio", 0.0),
            "fws": fw, "hw": hw,
            "mean_depth": depth if depth is not None else 0.0,
            "heterozygosity": het,
        }
        v = judge(metrics, thr)
        rep_dir = Path(config.output_dir) / sample / "08_report"
        rep_dir.mkdir(parents=True, exist_ok=True)
        paths = {
            "vaf_histogram": f"../05_vaf/{sample}.vaf_histogram.png",
            "genomescope": f"../../07_kmer/{sample}/02_genomescope/linear_plot.png",
        }
        (rep_dir / f"{sample}.report.md").write_text(build_sample_report(sample, metrics, v, paths))
        rows.append({
            "sample": sample, "verdict": v["verdict"], "confidence": v["confidence"],
            "vaf_mid_ratio": round(metrics["vaf_mid_ratio"], 4),
            "multiallelic_ratio": round(metrics["multiallelic_ratio"], 4),
            "fws": round(fw, 3) if fw is not None else "N/A",
            "mean_depth": round(metrics["mean_depth"], 1),
        })
        logger.info(f"{sample}: {v['verdict']} (置信|confidence {v['confidence']}) — {v['votes']}")
    return rows


def main():
    """主入口|main entry."""
    config = _argv_to_config()
    config.validate()
    log_file = str(Path(config.output_dir) / "99_logs" / "mixrace.log")
    logger = ModuleLogger(log_file=log_file).get_logger()
    logger.info(f"mixrace 启动|mixrace start (fastq_dir={config.fastq_dir}, step={config.step})")
    runner = CommandRunner(logger, dry_run=config.dry_run)
    ckpt_dir = str(Path(config.output_dir) / "00_pipeline_info" / "checkpoints")
    ckpt = CheckpointManager(ckpt_dir, logger)
    try:
        run_pipeline(config, runner, ckpt, logger)
    except Exception as e:
        logger.error(f"流程失败|pipeline failed: {e}")
        raise
    logger.info("mixrace 完成|mixrace done")


if __name__ == "__main__":
    main()
