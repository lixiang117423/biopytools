"""mixrace 主入口与编排(v0.3 GTX 后端)|mixrace main (GTX backend, three-branch verdict).

解析参数 → 按 --step 调度:01 QC+寄主剔除 / 02 GTX 比对+联合calling /
03 杂合评估+三分支判读 / 04 mapped-reads k-mer / 05 图+报告。
|argparse -> --step dispatch: QC+host depletion / GTX / het-eval verdicts /
mapped-read k-mer / figures+report. Mentor methodology v4.
"""
import argparse
import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace as dc_replace
from pathlib import Path

from .config import MixraceConfig
from .utils import ModuleLogger, CommandRunner, CheckpointManager, write_software_versions
from .samples import discover_samples
from .pipeline import run_index, run_qc, run_depth, read_cached_depth, run_kmer
from .host_filter import run_host_index, run_host_filter, pathogen_alignment_stats
from .gtx_backend import run_gtx, extract_mapped_fastq, count_mapped
from .het_eval import run_het_eval, write_tsv


def _worker_config(config):
    """样本级并行时每 worker 的线程配额(t/并行数)|per-worker thread quota."""
    n = getattr(config, "sample_parallel", 1)
    if n <= 1:
        return config
    return dc_replace(config, threads=max(1, config.threads // n))


def _parallel_map(config, runner, items, fn):
    """逐样本并行(默认1=串行);每 worker 独立 CommandRunner(避免 _child_proc 互踩)。
    |per-sample parallelism; each worker owns its CommandRunner."""
    n = getattr(config, "sample_parallel", 1)
    if n <= 1 or len(items) <= 1:
        return [fn(_worker_config(config), runner, it) for it in items]

    def _one(it):
        w_runner = CommandRunner(runner.logger, dry_run=runner.dry_run)
        return fn(_worker_config(config), w_runner, it)
    with ThreadPoolExecutor(max_workers=n) as ex:
        return list(ex.map(_one, items))


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
                                description="WGS混合小种检测(三分支判读)|WGS mixed-race detection")
    p.add_argument("-i", "--input", dest="fastq_dir", default=None)
    p.add_argument("--clean-fastq-dir", dest="clean_fastq_dir", default=None)
    p.add_argument("-g", "--genome", required=True)
    p.add_argument("-o", "--output-dir", dest="output_dir", default="mixrace_out")
    p.add_argument("--repeat-bed", dest="repeat_bed", default=None,
                   help="额外排除区域BED(与自动热点并集)|extra exclude BED (merged with hotspots)")
    p.add_argument("--host-genome", dest="host_genome", default=None,
                   help="寄主基因组FASTA(给则比对寄主并整对剔除寄主reads,报告寄主占比)"
                        "|host genome FASTA (deplete host reads, report host rate)")
    p.add_argument("--min-mapq", dest="min_mapq", type=int, default=20,
                   help="比对质量阈值:mapped reads提取+统计口径(0=不过滤)"
                        "|min MAPQ for mapped-read extraction + stats (0=off)")
    p.add_argument("-t", "--threads", type=int, default=12)
    p.add_argument("--sample-parallel", dest="sample_parallel", type=int, default=1,
                   help="样本级并行数(寄主剔除/mapped提取/reads统计同时跑N个样本,"
                        "每worker线程=threads/N;默认1串行)|per-sample parallelism")
    p.add_argument("-k", "--kmer-size", type=int, default=21)
    p.add_argument("-l", "--read-length", type=int, default=150)
    p.add_argument("--step", type=int, default=None)
    p.add_argument("--no-checkpoint", dest="enable_checkpoint", action="store_false")
    p.add_argument("--dry-run", dest="dry_run", action="store_true")
    # 三分支判读阈值|three-branch verdict thresholds
    p.add_argument("--pure-het-threshold", dest="pure_het_threshold", type=float, default=0.001,
                   help="总杂合率低于此值判纯菌(默认0.001=0.1%%)|pure threshold")
    p.add_argument("--partner-alt-rate", dest="partner_alt_min", type=float, default=0.8,
                   help="混合伴侣:ALT携带率阈值(默认0.8)|partner ALT-carrier threshold")
    p.add_argument("--partner-hom-rate", dest="partner_hom_min", type=float, default=0.5,
                   help="混合伴侣:伴侣纯合1/1占比阈值(默认0.5)|partner homozygous threshold")
    p.add_argument("--min-sites", dest="min_sites", type=int, default=1000,
                   help="最低有GT位点数,低于判uncertain(默认1000)|min called sites")
    p.add_argument("--window-size", dest="window_size", type=int, default=100000,
                   help="热点窗口大小bp(默认100kb)|hotspot window size")
    p.add_argument("--hotspot-fold", dest="hotspot_fold", type=float, default=2.0,
                   help="热点:窗口杂合率>该倍数×自身全基因组率(默认2)|hotspot fold")
    p.add_argument("--hotspot-min-median", dest="hotspot_min_median", type=float, default=0.10,
                   help="热点:窗口在候选中的中位杂合率下限(默认0.1)|hotspot min median rate")
    a = p.parse_args()
    return MixraceConfig(
        fastq_dir=a.fastq_dir or "", clean_fastq_dir=a.clean_fastq_dir,
        genome=a.genome, output_dir=a.output_dir, repeat_bed=a.repeat_bed,
        host_genome=a.host_genome, min_mapq=a.min_mapq,
        pure_het_threshold=a.pure_het_threshold, partner_alt_min=a.partner_alt_min,
        partner_hom_min=a.partner_hom_min, min_sites=a.min_sites,
        window_size=a.window_size, hotspot_fold=a.hotspot_fold,
        hotspot_min_median=a.hotspot_min_median,
        threads=a.threads, sample_parallel=a.sample_parallel,
        kmer_size=a.kmer_size, read_length=a.read_length,
        step=a.step, enable_checkpoint=a.enable_checkpoint, dry_run=a.dry_run)


def _read_verdict_table(config) -> list:
    """--step 4/5 重跑:读已判读表(实现在 het_eval)|reread verdict table (impl in het_eval)."""
    from .het_eval import read_verdict_table
    return read_verdict_table(config)


def _reads_accounting(config, runner, rows, bam_dir: str, genome_size: int):
    """step3b: 逐样本 reads 账本(host/mapping/污染)+深度+覆盖广度。
    |per-sample reads accounting (host/map/contamination) + depth + breadth."""
    eval_dir = Path(config.output_dir) / "04_het_eval"
    eval_dir.mkdir(parents=True, exist_ok=True)

    def _account_one(wcfg, w_runner, row):
        sample = row["sample"]
        bam = str(Path(bam_dir) / f"{sample}.bam")
        bam_exists = Path(bam).exists()
        # 计数表(供 pathogen_alignment_stats 读)|counts file for pathogen stats
        if bam_exists and not (eval_dir / f"{sample}.mapq_stats.tsv").exists():
            total, mapped = count_mapped(w_runner, wcfg, bam)
            if total is not None and mapped is not None:
                (eval_dir / f"{sample}.mapq_stats.tsv").write_text(
                    f"field\tvalue\ntotal_primary_reads\t{total}\n"
                    f"mapped_q_reads\t{mapped}\n", encoding="utf-8")
        stats_file = eval_dir / "alignment_qc" / f"{sample}.stats.txt"
        depth = read_cached_depth(stats_file, genome_size)
        if depth is None and bam_exists:
            depth = run_depth(w_runner, wcfg, sample, bam, genome_size)
        depth = depth if depth is not None else 0.0
        try:
            hs = pathogen_alignment_stats(wcfg, w_runner, sample, bam, genome_size,
                                          mean_depth=depth)
        except Exception as e:
            w_runner.logger.warning(f"{sample}: reads统计失败|reads stats failed: {e}")
            hs = None
        return sample, depth, hs

    stats_by_sample = {}
    for sample, depth, hs in _parallel_map(config, runner, rows, _account_one):
        stats_by_sample[sample] = hs
        row = next(r for r in rows if r["sample"] == sample)
        row["mean_depth"] = depth
        if hs:
            row["host_rate"] = hs.get("host_rate")
            row["pathogen_map_rate"] = hs.get("pathogen_map_rate")
            row["contamination_rate"] = hs.get("unassigned_rate")
            row["breadth_1x"] = hs.get("breadth_1x")
    return stats_by_sample


def _figures_and_report(config, runner, ckpt, logger, rows, genome_size):
    """step5: 全套图 + 判读汇总表 + 单样本报告 + HTML|figures + reports."""
    from .reporter import build_sample_report, build_summary_table, build_html_report
    rep_dir = Path(config.output_dir) / "07_report"
    rep_dir.mkdir(parents=True, exist_ok=True)
    summ_dir = Path(config.output_dir) / "summary"
    summ_dir.mkdir(parents=True, exist_ok=True)
    # 图(失败降级,不阻断报告)|figures (degrade gracefully)
    figures = {}
    try:
        from .figures import build_figures
        fig_dir = Path(config.output_dir) / "06_figures"
        payloads = {"rows": rows,
                    "het_eval_dir": Path(config.output_dir) / "04_het_eval"}
        paths = build_figures(config, logger, fig_dir, payloads)
        figures = {p.stem: str(p) for p in paths}
    except Exception as e:
        logger.warning(f"绘图失败(报告继续)|figure generation failed: {e}")
    # 汇总表|summary table
    tsv, html = build_summary_table(rows)
    (summ_dir / "verdict_summary.tsv").write_text(tsv, encoding="utf-8")
    (summ_dir / "verdict_summary.html").write_text(html, encoding="utf-8")
    # 单样本 md|per-sample markdown
    for row in rows:
        (rep_dir / f"{row['sample']}.report.md").write_text(
            build_sample_report(row["sample"], row, figures), encoding="utf-8")
    # 自包含 HTML|self-contained HTML
    # 判读口径文案按 config 阈值动态生成(防阈值可配后文案漂移)
    # |verdict note built from config thresholds (no drift)
    verdict_note = (f"判读口径: 杂合率<{config.pure_het_threshold*100:.1f}% 纯菌;"
                    f"强混合伴侣(ALT携带≥{config.partner_alt_min*100:.0f}%且纯合1/1≥"
                    f"{config.partner_hom_min*100:.0f}%)=混杂菌株;其余=优势菌株/参考差异型。"
                    f"建议列为实验操作指引(可保存/需再分离纯化)。")
    (summ_dir / "mixrace_report.html").write_text(
        build_html_report("根肿菌样本混杂评估报告", rows, figures, verdict_note),
        encoding="utf-8")
    logger.info(f"汇总表已写|summary written: {summ_dir / 'verdict_summary.tsv'}")


def run_pipeline(config, runner, ckpt, logger):
    """5 步编排|orchestrate 5 steps by --step."""
    step = config.step
    genome_size = _genome_size(config.genome)
    logger.info(f"基因组大小|genome size: {genome_size} bp")

    # 01 QC + 寄主剔除|QC + host depletion
    if step in (None, 1):
        run_index(config, runner, ckpt)
        clean_dir = run_qc(config, runner, ckpt)
    else:
        clean_dir = (Path(config.clean_fastq_dir) if config.clean_fastq_dir
                     else Path(config.output_dir) / "01_qc")

    samples = []
    if step in (None, 1, 2, 4):
        src = str(clean_dir) if Path(str(clean_dir)).is_dir() else config.fastq_dir
        samples = discover_samples(src) if src and os.path.isdir(src) else []
        logger.info(f"检测到 {len(samples)} 个样本|{len(samples)} samples detected")

    # 01b 寄主剔除(--host-genome)|host depletion
    host_failed = set()
    depleted_dir = clean_dir
    if config.host_genome:
        host_dir = Path(config.output_dir) / "02_host_filter"
        if step in (None, 1):
            host_idx = run_host_index(config, runner, ckpt)

            def _host_one(wcfg, w_runner, s):
                return run_host_filter(wcfg, w_runner, ckpt, s["sample"],
                                       s["r1"], s["r2"], str(host_idx))
            results = _parallel_map(config, runner, samples, _host_one)
            for s, res in zip(samples, results):
                if res:
                    s["r1"], s["r2"] = res["nohost_r1"], res["nohost_r2"]
                else:
                    host_failed.add(s["sample"])
            depleted_dir = host_dir   # 真实 run_host_filter 已建目录|dir created by host filter
        elif host_dir.is_dir():
            # step>=2 单独重跑:从 nohost 产物重新发现样本|rediscover from nohost outputs
            samples = discover_samples(str(host_dir))
            depleted_dir = host_dir
            logger.info(f"寄主剔除目录发现 {len(samples)} 个样本|{len(samples)} samples "
                        f"from 02_host_filter dir")
        else:
            logger.warning("启用了寄主剔除但 02_host_filter/ 不存在(先跑 --step 1),"
                           "下游沿用 clean fastq|02_host_filter/ missing, using clean fastq")

    # 02 GTX 比对+联合calling|GTX mapping + joint calling
    bam_dir = str(Path(config.output_dir) / "03_gtx" / "03_mapping" / "bam")
    vcf = str(Path(config.output_dir) / "03_gtx" / "04_joint_calling" /
              "gtx_joint_raw.vcf.gz")
    if step in (None, 2):
        bam_dir, vcf = run_gtx(config, runner, ckpt, str(depleted_dir))
        if vcf is None:
            logger.error("GTX 失败,后续步骤中止|GTX failed, downstream aborted")
            write_software_versions(
                config, logger,
                str(Path(config.output_dir) / "00_pipeline_info" / "software_versions.yml"))
            return

    # 03 杂合评估+判读+reads账本|het eval + verdicts + reads accounting
    if step in (None, 3):
        rows = run_het_eval(config, runner, ckpt, vcf)
        if rows:
            _reads_accounting(config, runner, rows, bam_dir, genome_size)
    else:
        rows = _read_verdict_table(config)
        if step == 5:
            _reads_accounting(config, runner, rows, bam_dir, genome_size)

    # 04 mapped reads 提取 + k-mer|mapped reads + smudgescope
    if step in (None, 4):
        mapped_dir = Path(config.output_dir) / "05_kmer" / "mapped_fastq"
        names = [r["sample"] for r in rows] or [s["sample"] for s in samples]
        todo = [n for n in names if n not in host_failed
                and (runner.dry_run or (Path(bam_dir) / f"{n}.bam").exists())]
        for n in set(names) - set(todo):
            if n not in host_failed:
                logger.warning(f"{n}: BAM 缺失,跳过 k-mer 输入|BAM missing, k-mer skipped")

        def _mapped_one(wcfg, w_runner, name):
            return extract_mapped_fastq(wcfg, w_runner, ckpt, name,
                                        str(Path(bam_dir) / f"{name}.bam"))
        _parallel_map(config, runner, todo, _mapped_one)
        run_kmer(config, runner, ckpt, str(mapped_dir))

    # 05 图+报告|figures + report
    if step in (None, 5) and rows:
        _figures_and_report(config, runner, ckpt, logger, rows, genome_size)

    write_software_versions(
        config, logger,
        str(Path(config.output_dir) / "00_pipeline_info" / "software_versions.yml"))


def main():
    config = _argv_to_config()
    config.validate()
    log_file = str(Path(config.output_dir) / "99_logs" / "mixrace.log")
    logger = ModuleLogger(log_file=log_file).get_logger()
    logger.info(f"mixrace 启动|mixrace start (step={config.step}, "
                f"host_genome={'yes' if config.host_genome else 'no'})")
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
