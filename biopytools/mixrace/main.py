"""mixrace 主入口与编排(单倍体 freebayes 后端)|mixrace main (haploid freebayes backend).

解析参数 → 按 --step 调度:01 索引+QC / 02 比对+markdup / 03 freebayes -p 1 / 04 过滤 /
05 AFS(杂合率+形态+优势株占比) / 06 k-mer(smudgescope) / 07 判读+报告。逐样本。
|argparse -> --step dispatch over 7 per-sample steps (bwa-mem2+markdup, freebayes -p 1,
AFS analysis, smudgescope, verdict). Advisor methodology v2.
"""
import argparse
import os
import sys
from pathlib import Path

from .config import MixraceConfig
from .utils import ModuleLogger, CommandRunner, CheckpointManager, write_software_versions
from .samples import discover_samples
from .pipeline import (run_index, run_qc, run_align, run_call_freebayes, run_filter,
                       run_depth, run_kmer, parse_genomescope_model)
from .vaf_analysis import parse_freebayes, compute_afs
from .verdict import judge, calibrate_thresholds
from .reporter import (generate_vaf_histogram_r, build_sample_report, build_summary_table,
                       build_html_report)


def _genome_size(fa_path: str) -> int:
    """纯 Python 统计基因组大小(非 header 行的非空白字符数;支持 .gz)|
    genome size via base count (gzip-aware: plain open on .gz counts binary junk)."""
    import gzip
    n = 0
    try:
        opener = gzip.open if fa_path.endswith(".gz") else open
        with opener(fa_path, "rt") as fh:
            for line in fh:
                if not line.startswith(">"):
                    n += len(line.strip())
    except OSError:
        return 0
    return n


def _argv_to_config() -> MixraceConfig:
    p = argparse.ArgumentParser(prog="mixrace",
                                description="WGS混合小种检测(单倍体)|WGS mixed-race detection (haploid)")
    p.add_argument("-i", "--input", dest="fastq_dir", default=None)
    p.add_argument("--clean-fastq-dir", dest="clean_fastq_dir", default=None)
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
    p.add_argument("--min-coverage", type=int, default=30, help="freebayes --min-coverage(默认30)")
    p.add_argument("--min-alt-fraction", type=float, default=0.02,
                   help="freebayes --min-alternate-fraction(默认0.02,保低频等位)")
    p.add_argument("--pure-samples", dest="pure_samples", default=None,
                   help="已知纯样品(逗号分隔,校准het阈值)|known-pure samples (calibrate)")
    a = p.parse_args()
    pure = a.pure_samples.split(",") if a.pure_samples else None
    return MixraceConfig(
        fastq_dir=a.fastq_dir or "", clean_fastq_dir=a.clean_fastq_dir,
        genome=a.genome, output_dir=a.output_dir, repeat_bed=a.repeat_bed,
        threads=a.threads, kmer_size=a.kmer_size, read_length=a.read_length,
        step=a.step, enable_checkpoint=a.enable_checkpoint, dry_run=a.dry_run,
        min_qual=a.min_qual, min_dp=a.min_dp, min_alt_reads=a.min_alt_reads,
        freebayes_min_coverage=a.min_coverage,
        freebayes_min_alternate_fraction=a.min_alt_fraction, pure_samples=pure)


def _sample_afs(config, runner, sample: str, filt_vcf: str, genome_size: int):
    """从 freebayes VCF 查 AO/RO → 杂合率/AFS形态/优势株占比;写 vaf.tsv。
    |query AO/RO from freebayes VCF -> het_rate/AFS/dominant; write vaf.tsv."""
    ok, qtxt, _ = runner.run_conda(
        config.bcftools_path,
        ["query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%RO]\t[%AO]\n", str(filt_vcf)],
        f"query AO/RO|query AO/RO {sample}")
    recs = parse_freebayes(qtxt) if ok else []
    afs = compute_afs(recs, config.min_alt_reads, genome_size)
    # 写 vaf.tsv(每位点 RO/AO/VAF)|write per-site RO/AO/VAF
    vaf_dir = Path(config.output_dir) / "05_vaf"
    vaf_dir.mkdir(parents=True, exist_ok=True)
    vaf_tsv = vaf_dir / f"{sample}.vaf.tsv"
    with open(vaf_tsv, "w") as fh:
        fh.write("chrom\tpos\tref\talts\tro\taos\tvafs\n")
        for r in recs:
            tot = (r["ro"] + sum(r["aos"])) or 1
            fh.write("\t".join([
                r["chrom"], str(r["pos"]), r["ref"], ",".join(r["alts"]),
                str(r["ro"]), ",".join(str(x) for x in r["aos"]),
                ",".join(f"{a/tot:.4f}" for a in r["aos"])]) + "\n")
    return afs


def run_pipeline(config, runner, ckpt, logger):
    """7 步编排(单倍体,逐样本)|orchestrate 7 per-sample steps by --step."""
    step = config.step
    thr = {"het_pure": config.het_pure, "het_suspicious": config.het_suspicious,
           "het_impure": config.het_impure, "min_depth": config.min_depth}
    genome_size = _genome_size(config.genome)
    logger.info(f"基因组大小|genome size: {genome_size} bp")

    # 01 索引 + QC|index + QC
    if step in (None, 1):
        idx_dir = run_index(config, runner, ckpt)
        clean_dir = run_qc(config, runner, ckpt)
    else:
        idx_dir = Path(config.output_dir) / "00_pipeline_info" / "index"
        clean_dir = (Path(config.clean_fastq_dir) if config.clean_fastq_dir
                     else Path(config.output_dir) / "01_qc")

    samples = []
    if step in (None, 2, 3, 4, 5, 7):
        src = str(clean_dir) if Path(str(clean_dir)).is_dir() else config.fastq_dir
        samples = discover_samples(src) if src and os.path.isdir(src) else []
        logger.info(f"检测到 {len(samples)} 个样本|{len(samples)} samples detected")

    fa = str(Path(str(idx_dir)) / os.path.basename(config.genome))
    afs = {}
    # 02-05 逐样本:比对→freebayes→过滤→AFS|per-sample align→freebayes→filter→AFS
    for s in samples:
        sample = s["sample"]
        if step in (None, 2):
            bam = run_align(config, runner, ckpt, sample, s["r1"], s["r2"], str(idx_dir))
        else:
            bam = Path(config.output_dir) / "02_alignment" / f"{sample}.sorted.markdup.bam"
        if step in (None, 3):
            vcf = run_call_freebayes(config, runner, ckpt, sample, str(bam), fa)
        else:
            vcf = Path(config.output_dir) / "03_variants" / f"{sample}.raw.vcf.gz"
        if step in (None, 4):
            filt = run_filter(config, runner, ckpt, sample, str(vcf))
        else:
            filt = Path(config.output_dir) / "04_filtered" / f"{sample}.filtered.vcf.gz"
        if step in (None, 5):
            afs[sample] = _sample_afs(config, runner, sample, str(filt), genome_size)
            vaf_tsv = Path(config.output_dir) / "05_vaf" / f"{sample}.vaf.tsv"
            png = Path(config.output_dir) / "05_vaf" / f"{sample}.vaf_histogram.png"
            script = generate_vaf_histogram_r(str(vaf_tsv), str(png), config.rscript_path)
            runner.run_conda(config.rscript_path, [script], f"AFS直方图|AFS histogram {sample}")

    # 06 k-mer 谱(smudgescope,读 clean fastq)|k-mer spectrum via smudgescope
    if step in (None, 6):
        run_kmer(config, runner, ckpt, str(clean_dir))

    # 07 判读 + 报告|verdict + report
    if step in (None, 7):
        rows = _verdict_and_report(config, runner, logger, samples, afs, thr, genome_size)
        summ_dir = Path(config.output_dir) / "summary"
        summ_dir.mkdir(parents=True, exist_ok=True)
        tsv, html = build_summary_table(rows)
        (summ_dir / "verdict_summary.tsv").write_text(tsv)
        (summ_dir / "verdict_summary.html").write_text(html)
        logger.info(f"汇总表已写|summary written: {summ_dir / 'verdict_summary.tsv'}")

    write_software_versions(
        config, logger,
        str(Path(config.output_dir) / "00_pipeline_info" / "software_versions.yml"))


def _read_heterozygosity(config, sample: str):
    model = Path(config.output_dir) / "06_kmer" / sample / "02_genomescope" / "model.txt"
    if model.exists():
        return parse_genomescope_model(model.read_text()).get("heterozygosity")
    return None


def _verdict_and_report(config, runner, logger, samples, afs, thr, genome_size):
    """step07: 逐样品判读 + 报告(支持单独重跑:缺 AFS 则从已过滤 VCF 重算)。|step07 verdict+report."""
    filt_base = Path(config.output_dir)
    for sample in [s["sample"] for s in samples]:
        if sample not in afs:
            filt = filt_base / "04_filtered" / f"{sample}.filtered.vcf.gz"
            if filt.exists():
                afs[sample] = _sample_afs(config, runner, sample, str(filt), genome_size)

    if config.pure_samples:
        pure_rows = [afs[p] for p in config.pure_samples if p in afs]
        if len(pure_rows) >= 2:
            cal = calibrate_thresholds(pure_rows)
            if cal.get("het_pure") is not None:
                thr["het_pure"] = cal["het_pure"]
            logger.info(f"het 阈值已用 {len(pure_rows)} 个纯样品校准|het threshold calibrated")

    rows = []
    samples_data = []
    for s in samples:
        sample = s["sample"]
        a = afs.get(sample, {})
        # 数据缺失(VCF缺/query失败/0变异)≠真纯,交 judge() no_data 走 uncertain|
        # Missing data (no VCF/query failed/0 variants) != pure; judge() no_data -> uncertain.
        no_data = (not a) or a.get("total_variant", 0) == 0
        if no_data:
            logger.warning(f"{sample}: 未获取到变异数据,判读 uncertain(先跑 --step 5)"
                           f"|No variant data, verdict uncertain (run --step 5 first)")
        bam = filt_base / "02_alignment" / f"{sample}.sorted.markdup.bam"
        depth = run_depth(runner, config, sample, str(bam), genome_size) if bam.exists() else None
        het = _read_heterozygosity(config, sample)
        metrics = {
            "het_rate": a.get("het_rate", 0.0), "het_sites": a.get("het_sites", 0),
            "afs_shape": a.get("afs_shape", "monoclonal"),
            "dominant_proportion": a.get("dominant_proportion"),
            "maf": a.get("maf"),
            "mean_depth": depth if depth is not None else 0.0,
            "heterozygosity": het,
        }
        v = judge({**metrics, "no_data": no_data}, thr)
        # 纯样品不显示优势株占比(杂合位点太少,占比无统计意义)|hide dominant for pure samples
        show_dom = None if v["verdict"] == "single_genotype" else metrics["dominant_proportion"]
        metrics_disp = {**metrics, "dominant_proportion": show_dom}
        rep_dir = Path(config.output_dir) / "07_report"
        rep_dir.mkdir(parents=True, exist_ok=True)
        paths = {"afs_histogram": f"../05_vaf/{sample}.vaf_histogram.png",
                 "genomescope": f"../06_kmer/{sample}/02_genomescope/linear_plot.png"}
        (rep_dir / f"{sample}.report.md").write_text(build_sample_report(sample, metrics_disp, v, paths))
        rows.append({
            "sample": sample, "verdict": v["verdict"], "confidence": v["confidence"],
            "het_rate": f"{metrics['het_rate']*100:.4f}%",
            "afs_shape": metrics["afs_shape"],
            "dominant_proportion": f"{show_dom*100:.1f}%" if show_dom is not None else "—",
            "mean_depth": round(metrics["mean_depth"], 1),
        })
        # 合并 HTML 报告用的完整数据(含图片绝对路径)|full data + image paths for merged HTML
        samples_data.append({
            "sample": sample, "verdict": v["verdict"], "confidence": v["confidence"],
            "rationale": v["rationale"],
            "het_rate": metrics["het_rate"], "afs_shape": metrics["afs_shape"],
            "dominant_proportion": show_dom, "mean_depth": metrics["mean_depth"],
            "metrics": metrics_disp,
            "images": {
                "afs": str(filt_base / "05_vaf" / f"{sample}.vaf_histogram.png"),
                "genomescope": str(filt_base / "06_kmer" / sample / "02_genomescope" / "linear_plot.png"),
                "smudgeplot": str(filt_base / "06_kmer" / sample / "03_smudgeplot" / f"{sample}_smudgeplot.png"),
            },
        })
        logger.info(f"{sample}: {v['verdict']} (置信|confidence {v['confidence']}) "
                    f"het={metrics['het_rate']*100:.4f}% shape={metrics['afs_shape']}")
    # 合并 HTML 报告(单文件,图片内嵌)|merged self-contained HTML report
    html_path = Path(config.output_dir) / "summary" / "mixrace_report.html"
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(build_html_report("根肿菌混合小种检测报告", samples_data), encoding="utf-8")
    logger.info(f"HTML 报告已写|HTML report written: {html_path}")
    return rows


def main():
    config = _argv_to_config()
    config.validate()
    log_file = str(Path(config.output_dir) / "99_logs" / "mixrace.log")
    logger = ModuleLogger(log_file=log_file).get_logger()
    logger.info(f"mixrace 启动|mixrace start (step={config.step}, clean_fastq_dir={config.clean_fastq_dir})")
    runner = CommandRunner(logger, dry_run=config.dry_run)
    ckpt = CheckpointManager(str(Path(config.output_dir) / "00_pipeline_info" / "checkpoints"), logger)
    try:
        run_pipeline(config, runner, ckpt, logger)
    except Exception as e:
        logger.error(f"流程失败|pipeline failed: {e}")
        raise
    logger.info("mixrace 完成|mixrace done")


if __name__ == "__main__":
    main()
