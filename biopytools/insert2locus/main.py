"""insert2locus命令行入口|insert2locus CLI entry point

编排:01比对→02junction钓取→03迭代步移→04完整locus→05覆盖验证→06HTML报告|
Orchestration: 01 mapping -> 02 junction fishing -> 03 walking -> 04 locus
-> 05 verification -> 06 HTML report
"""

import argparse
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from . import data_processing, locus_builder, report, verifier, walking
from .config import Insert2locusConfig
from .data_processing import extract_flank_bait
from .locus_builder import CompleteLocus, annotate_junctions, \
    classify_anchor, find_complete_locus, insert_regions_from_alignment, \
    max_matched_by_query, write_junction_outputs
from .utils import CommandRunner, Insert2locusLogger, build_conda_command, \
    format_number
from .verifier import classify, count_junction_reads, per_record_depth, \
    segment_coverage
from .walking import WalkingRunner, contig_stats, fish_mates


def parse_arguments(argv=None):
    """解析参数|Parse arguments"""
    parser = argparse.ArgumentParser(
        prog="biopytools insert2locus",
        description="转基因插入位点提取:soft-clip钓取+迭代步移+完整locus重构+覆盖验证"
                    "|Transgenic insertion locus extraction",
        epilog="示例|Examples: biopytools insert2locus -i fq_dir/ -f insert.fasta -o output/")
    parser.add_argument("-i", "--input", required=True,
                        help="fastq目录或R1文件|fastq dir or R1 file")
    parser.add_argument("-f", "--insert-fasta", required=True,
                        help="插入序列fasta(载体+片段)|Insert fasta (vector+fragment)")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="输出目录|Output dir")
    parser.add_argument("-t", "--threads", type=int, default=12,
                        help="线程数|Threads (default 12)")
    parser.add_argument("--sort-mem", default="2G",
                        help="samtools sort单线程内存|per-thread sort memory")
    parser.add_argument("--read1-suffix", default=None,
                        help="R1后缀(默认自动检测)|R1 suffix (auto-detect by default)")
    parser.add_argument("--max-rounds", type=int, default=30,
                        help="步移最大轮数|Max walking rounds")
    parser.add_argument("--min-softclip", type=int, default=25,
                        help="诱饵softclip最短长度|Min softclip for bait")
    parser.add_argument("--min-unmapped", type=int, default=400,
                        help="未比对contig作诱饵最短长度|Min unmapped length for bait")
    parser.add_argument("--min-growth", type=int, default=50,
                        help="收敛判定增量阈值|Growth threshold for convergence")
    parser.add_argument("--mapq-min", type=int, default=1,
                        help="招募最低MAPQ|Min MAPQ for recruitment")
    parser.add_argument("--repeat-cap", type=int, default=10000,
                        help="单轮新招reads上限(撞重复区)|Per-round new-read cap")
    parser.add_argument("--junction-flank", type=int, default=50,
                        help="边界报告最短侧翼|Min flank for border report")
    parser.add_argument("--tdna-fasta", default=None,
                        help="单独插入序列fasta(区分insert与载体骨架;不给则-f整体当insert)"
                             "|Standalone T-DNA fasta (separates insert from backbone)")
    parser.add_argument("--target-flank", type=int, default=None,
                        help="LB/RB目标侧翼长度(默认不限,尽可能走远到自然收敛;"
                             "到达后小增量即收敛)"
                             "|Target flank length (default: walk as far as possible)")
    parser.add_argument("--force", action="store_true",
                        help="忽略断点全部重跑|Ignore checkpoints and rerun")
    parser.add_argument("--log-level", default="INFO", help="日志级别|Log level")
    return parser.parse_args(argv)


def run_stage_with_checkpoint(name: str, markers: List, force: bool,
                              func, logger=None):
    """步级断点续传|Stage-level checkpoint resume"""
    if markers and not force and all(Path(m).exists() for m in markers):
        if logger:
            logger.info(f"跳过已完成步骤|Skipping completed step: {name}")
        return True
    if logger:
        logger.info(f"开始步骤|Starting step: {name}")
    return func()


# ---------- 通用小工具|Common helpers ----------

def read_fasta_dict(fasta: Path) -> Dict[str, str]:
    """读fasta成{name:seq}|Read fasta into {name:seq}"""
    seqs: Dict[str, str] = {}
    name = None
    chunks: List[str] = []
    for line in Path(fasta).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def read_insert_record_lens(insert_fasta: Path, work_dir: Path,
                            samtools_path: Optional[str] = None) -> Dict[str, int]:
    """读insert各record长度(.fai优先,缺失时用samtools生成)|
    Record lengths from .fai (generate via samtools when missing)"""
    fai = Path(str(insert_fasta) + ".fai")
    if not fai.exists() and samtools_path:
        CommandRunner(_silent_logger(), work_dir).run(
            [samtools_path, "faidx", str(insert_fasta)],
            description="faidx索引insert|faidx insert")
    if fai.exists():
        lens = {}
        for line in fai.read_text().splitlines():
            f = line.split("\t")
            if len(f) >= 2:
                lens[f[0]] = int(f[1])
        return lens
    # 无samtools时容错:纯Python读fasta|Fallback: pure-Python fasta parse
    return {name: len(seq) for name, seq in read_fasta_dict(insert_fasta).items()}


def _silent_logger():
    """无handler的logger(faidx容错路径用,消息自然丢弃)|
    Handler-less logger for the faidx fallback path (messages go nowhere)"""
    import logging
    return logging.getLogger("insert2locus")


def probe_tool_versions(cfg, runner) -> Dict[str, str]:
    """探测工具版本|Probe tool versions"""
    versions = {}
    probes = [
        ("bwa", build_conda_command(cfg.bwa_path, [])),
        ("samtools", build_conda_command(cfg.samtools_path, ["--version"])),
        ("seqkit", build_conda_command(cfg.seqkit_path, ["version"])),
        ("spades", build_conda_command(cfg.spades_path, ["--version"])),
    ]
    for name, cmd in probes:
        out = runner.run_capture_stderr(cmd, description=f"版本探测|Version: {name}")
        versions[name] = (out or "").strip().splitlines()[0][:120] if out else "unknown"
    return versions


# ---------- 阶段函数|Stage functions ----------

def run_mapping_stage(cfg, logger, runner, sample: str, r1: Path, r2: Path,
                      sdir: Path) -> Dict[str, Path]:
    """01:比对到insert伪基因组|01: map to insert pseudo-genome"""
    mdir = sdir / "01_mapping"
    bam = mdir / f"{sample}.vs_insert.sorted.bam"
    flagstat = mdir / f"{sample}.vs_insert.flagstat.txt"
    cov_tsv = mdir / f"{sample}.insert_coverage.tsv"

    if not Path(str(cfg.insert_fasta) + ".bwt").exists():
        if not runner.run([cfg.bwa_path, "index", str(cfg.insert_fasta)],
                          description="bwa索引insert|bwa index insert"):
            raise RuntimeError("bwa index失败|bwa index failed")
    rg = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
    ok = runner.run_pipeline(
        [[cfg.bwa_path, "mem", "-t", str(cfg.threads), "-R", rg,
          str(cfg.insert_fasta), str(r1), str(r2)],
         [cfg.samtools_path, "sort", "-@", str(cfg.threads),
          "-m", cfg.sort_mem, "-o", str(bam), "-"]],
        description="WGS比对到insert|Map WGS to insert")
    if not ok:
        raise RuntimeError("比对失败|Mapping failed")
    if not runner.run([cfg.samtools_path, "index", str(bam)],
                      description="索引bam|Index bam"):
        raise RuntimeError("bam索引失败|Index failed")
    if not runner.run([cfg.samtools_path, "flagstat", str(bam)],
                      description="flagstat", stdout_file=flagstat):
        raise RuntimeError("flagstat失败|flagstat failed")

    # 逐record覆盖(backbone线索)|Per-record coverage (backbone signal)
    depth_out = runner.run_capture(
        [cfg.samtools_path, "depth", "-a", str(bam)],
        description="insert深度统计|Insert depth")
    rows = []
    if depth_out:
        for line in depth_out.splitlines():
            f = line.split("\t")
            if len(f) >= 3:
                rows.append((f[0], int(f[1]), int(f[2])))
    total_lens = read_insert_record_lens(Path(cfg.insert_fasta), sdir, cfg.samtools_path)
    stats = per_record_depth(rows, total_lens)
    with open(cov_tsv, "w") as fh:
        fh.write("record\tlength\tmean_depth\tcovered_frac\n")
        for rec, s in stats.items():
            fh.write(f"{rec}\t{s['length']}\t{s['mean']}\t{s['covered_frac']}\n")

    # tdna分区覆盖(insert区/骨架区,backbone信号)|Region coverage via tdna
    record_depth_for_report: Dict[str, dict] = dict(stats)
    anchor_label = None
    insert_regions = _compute_insert_regions(cfg, runner, sdir)
    if insert_regions is not None:
        region_stats = _region_stats(cfg, runner, bam, sdir, insert_regions,
                                     total_lens, rows)
        with open(cov_tsv, "a") as fh:
            fh.write("#region\tlength\tmean_depth\tcovered_frac\n")
            for label, s in region_stats.items():
                fh.write(f"{label}\t{s['length']}\t{s['mean']}\t"
                         f"{s['covered_frac']}\n")
        record_depth_for_report = region_stats
        anchor_label = "insert区|insert"
        logger.info(
            f"insert区{_merge_spans(insert_regions)}"
            f"(共{sum(e - s + 1 for s, e in _merge_spans(insert_regions))}bp),"
            f"骨架区覆盖{region_stats['骨架区|backbone']['covered_frac'] * 100:.1f}%"
            f"|Insert region located, backbone coverage computed")

    logger.info(
        f"01完成:比对统计见flagstat;逐record覆盖{format_number(len(rows))}行|"
        f"Stage 01 done, per-record coverage written")
    return {"bam": bam, "flagstat": flagstat, "coverage_tsv": cov_tsv,
            "record_depth": record_depth_for_report,
            "anchor_label": anchor_label, "insert_regions": insert_regions}


def _region_stats(cfg, runner, bam, sdir: Path, insert_regions, total_lens,
                  rows=None):
    """insert区/骨架区分区覆盖|Insert-region vs backbone region coverage"""
    from .verifier import region_depth
    if rows is None:
        # 续跑重算:重新取深度行|Resume path: recompute depth rows
        depth_out = runner.run_capture(
            [cfg.samtools_path, "depth", "-a", str(bam)],
            description="insert深度统计(续跑)|Insert depth (resume)")
        rows = []
        if depth_out:
            for line in depth_out.splitlines():
                f = line.split("\t")
                if len(f) >= 3:
                    rows.append((f[0], int(f[1]), int(f[2])))
    rec_len = sum(total_lens.values())
    merged = _merge_spans(insert_regions)
    complement = _complement_spans(merged, 1, rec_len)
    return region_depth(rows, {"insert区|insert": merged,
                               "骨架区|backbone": complement})


def _compute_insert_regions(cfg, runner, work_dir: Path):
    """tdna比对构建→insert区段(无tdna返回None)|tdna vs construct -> insert regions"""
    if not cfg.tdna_fasta:
        return None
    out = runner.run_pipeline_capture(
        [[cfg.bwa_path, "mem", "-t", str(cfg.threads), str(cfg.insert_fasta),
          str(cfg.tdna_fasta)],
         [cfg.samtools_path, "view", "-"]],
        description="tdna定位到构建|Locate tdna on construct")
    if out is None:
        raise RuntimeError("tdna定位失败|tdna localization failed")
    regions = insert_regions_from_alignment(out.splitlines())
    if not regions:
        raise RuntimeError(
            "tdna与构建无匹配,请检查-t与-f是否对应|tdna does not match construct")
    return regions


def _merge_spans(spans):
    """合并区段|Merge spans"""
    spans = sorted(spans)
    merged = [list(spans[0])]
    for start, end in spans[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def _complement_spans(spans, start: int, end: int):
    """区段补集|Complement of spans within [start,end]"""
    out = []
    cursor = start
    for s, e in spans:
        if s > cursor:
            out.append((cursor, s - 1))
        cursor = max(cursor, e + 1)
    if cursor <= end:
        out.append((cursor, end))
    return out


def run_junction_stage(cfg, logger, runner, sample: str,
                       bam: Path, r1: Path, r2: Path, sdir: Path) -> Dict[str, Path]:
    """02:soft-clip与mate-unmapped钓取|02: soft-clip & mate-unmapped fishing"""
    jdir = sdir / "02_junction_reads"
    sc_fq = jdir / f"{sample}.softclip.fastq"
    pe1 = jdir / f"{sample}.flank_candidates_R1.fastq"
    pe2 = jdir / f"{sample}.flank_candidates_R2.fastq"
    mate_bam = jdir / f"{sample}.mate_unmapped.bam"

    # soft-clip reads:CIGAR含S|Soft-clip reads: CIGAR with S
    sam_out = runner.run_capture(
        [cfg.samtools_path, "view", "-h", "-F", "4", "-q", str(cfg.mapq_min),
         str(bam)],
        description="提取softclip候选|Extract soft-clip candidates")
    if sam_out is None:
        raise RuntimeError("读取bam失败|Read bam failed")
    sc_sam = jdir / f"{sample}.softclip_reads.sam"
    with open(sc_sam, "w") as fh:
        for line in sam_out.splitlines():
            if line.startswith("@") or "S" in line.split("\t")[5]:
                fh.write(line + "\n")
    ok = runner.run_pipeline(
        [[cfg.samtools_path, "sort", "-@", str(cfg.threads), "-o",
          str(jdir / f"{sample}.softclip.sorted.bam"), str(sc_sam)]],
        description="softclip排序|Sort softclip")
    if not ok:
        raise RuntimeError("softclip排序失败|Sort softclip failed")
    fq_out = runner.run_capture(
        [cfg.samtools_path, "bam2fq",
         str(jdir / f"{sample}.softclip.sorted.bam")],
        description="softclip转fastq|softclip to fastq")
    if fq_out is None:
        raise RuntimeError("bam2fq失败|bam2fq failed")
    sc_fq.write_text(fq_out)

    # mate-unmapped reads|Mate-unmapped reads
    if not runner.run(
            [cfg.samtools_path, "view", "-b", "-f", "8", "-F", "260",
             str(bam)], description="mate未比对reads|Mate-unmapped reads",
            stdout_file=mate_bam):
        raise RuntimeError("mate_unmapped提取失败|Extract mate_unmapped failed")
    names_out = runner.run_capture(
        [cfg.samtools_path, "view", str(mate_bam)],
        description="读mate_unmapped|Read mate_unmapped")
    names = sorted({ln.split("\t")[0] for ln in (names_out or "").splitlines()
                    if ln and not ln.startswith("@")})
    logger.info(f"softclip/mate-unmapped钓取: {format_number(len(names))}条read名|"
                f"Junction fishing: {format_number(len(names))} read names")
    for mate, fastq, out in [(1, r1, pe1), (2, r2, pe2)]:
        if not fish_mates(cfg.seqkit_path, runner, names, mate, fastq,
                          jdir / f"{sample}.pat_mate{mate}.txt", out,
                          append=False):
            raise RuntimeError(f"钓mate/{mate}失败|Fish mate/{mate} failed")
    return {"softclip_fastq": sc_fq, "pe1": pe1, "pe2": pe2,
            "mate_unmapped_bam": mate_bam}


def run_locus_stage(cfg, logger, runner, sample: str, bam: Path, r1: Path,
                    r2: Path, sc_fq: Path, walk_dir: Path, sdir: Path,
                    insert_regions=None) -> dict:
    """04:招募全部覆盖reads重组装,判完整locus|04: recruit all, rebuild, find locus"""
    ldir = sdir / "04_locus"
    pool1 = ldir / f"{sample}.locus_pool_R1.fastq"
    pool2 = ldir / f"{sample}.locus_pool_R2.fastq"
    contigs = ldir / f"{sample}.contigs.fasta"
    if insert_regions is None:
        insert_regions = _compute_insert_regions(cfg, runner, sdir)

    # 招募全部比对上insert的reads名(含mate)|All mapped read names (mates included)
    # ⚠️ 不加-q过滤:MAPQ0多比对reads是拼通构建内部重复区必需的
    #   (人工分析用全部mapped 3911条才拼出贯穿contig,-q1只剩2025拼不出)|
    # No -q filter: MAPQ-0 multi-mappers are required to bridge construct
    # repeats (manual analysis used ALL 3911 mapped reads)
    # ⚠️ 不能加尾部"-":bam有.bai时samtools把后续参数当region解析→空输出|
    # No trailing "-": with an indexed bam samtools parses it as a region
    hits_out = runner.run_pipeline_capture(
        [[cfg.samtools_path, "view", "-F", "4", "-F", "260", str(bam)],
         ["cut", "-f1"],
         ["sort", "-u"]],
        description="招募insert覆盖reads|Recruit insert-covering reads")
    if hits_out is None:
        raise RuntimeError("招募失败|Recruit failed")
    names = [ln for ln in hits_out.splitlines() if ln]
    if not names:
        raise RuntimeError(
            "招募到0条reads:比对到构建的reads为空,请检查01产物|"
            "Recruited 0 reads: nothing maps to the construct, check stage 01")
    logger.info(f"招募全部覆盖reads(含MAPQ0): {format_number(len(names))}条|"
                f"Recruited {format_number(len(names))} read names (incl. MAPQ0)")
    for mate, fastq, out in [(1, r1, pool1), (2, r2, pool2)]:
        if not fish_mates(cfg.seqkit_path, runner, names, mate, fastq,
                          ldir / f"{sample}.pat_locus_mate{mate}.txt", out,
                          append=False):
            raise RuntimeError(f"钓locus mate/{mate}失败|Fish locus mate/{mate} failed")
    # 并入步移招募池|Merge walking recruited pools
    for pool, recruited in [(pool1, walk_dir / "recruited_R1.fastq"),
                            (pool2, walk_dir / "recruited_R2.fastq")]:
        if recruited.exists():
            with open(pool, "a") as fh:
                fh.write(recruited.read_text())

    # 终装:tdna作trusted-contigs骨架,穿过insert内部串联重复
    # (无tdna则纯de novo;实测rcrp内部重复会把contig劈成两半)|
    # Final assembly: tdna as trusted-contigs backbone to bridge internal
    # tandem repeats (without it the contig splits at the repeat)
    walker = WalkingRunner(cfg, logger, runner)
    trusted = Path(cfg.tdna_fasta) if cfg.tdna_fasta else None
    if not walker.run_spades(sc_fq, pool1, pool2, contigs,
                             trusted_contigs=trusted):
        raise RuntimeError("最终组装失败|Final assembly failed")
    contig_bam = ldir / f"{sample}.contigs_vs_insert.sorted.bam"
    if not walker.align_to_insert(contigs, contig_bam):
        raise RuntimeError("contigs比对失败|Contig alignment failed")
    sam_out = runner.run_capture(
        [cfg.samtools_path, "view", str(contig_bam)],
        description="读contigs比对|Read contig alignment")
    sam_lines = sam_out.splitlines() if sam_out else []
    contig_seqs = read_fasta_dict(contigs)

    # 完整locus锚定目标:有tdna锚tdna(区分骨架),否则锚构建|
    # Locus anchor: tdna when provided (distinguishes backbone), else construct
    if cfg.tdna_fasta:
        tdna_bam = ldir / f"{sample}.contigs_vs_tdna.sorted.bam"
        if not Path(str(cfg.tdna_fasta) + ".bwt").exists():
            if not runner.run([cfg.bwa_path, "index", str(cfg.tdna_fasta)],
                              description="bwa索引tdna|bwa index tdna"):
                raise RuntimeError("tdna index失败|tdna index failed")
        if not walker.align_to_reference(cfg.tdna_fasta, contigs, tdna_bam):
            raise RuntimeError("contigs比对tdna失败|Contigs vs tdna failed")
        locus_sam_out = runner.run_capture(
            [cfg.samtools_path, "view", str(tdna_bam)],
            description="读contigs_vs_tdna|Read contigs vs tdna")
        locus_sam = locus_sam_out.splitlines() if locus_sam_out else []
        anchor_lens = read_insert_record_lens(Path(cfg.tdna_fasta), sdir,
                                              cfg.samtools_path)
    else:
        locus_sam = sam_lines
        anchor_lens = read_insert_record_lens(Path(cfg.insert_fasta), sdir,
                                              cfg.samtools_path)
    locus = find_complete_locus(locus_sam, contig_seqs, anchor_lens)
    junctions = annotate_junctions(sam_lines, contig_seqs,
                                   cfg.junction_flank, cfg.min_unmapped,
                                   insert_regions=insert_regions)
    total, longest, n = contig_stats(contigs)
    logger.info(
        f"04完成:contigs={n} longest={longest}bp;完整locus="
        f"{'是|yes' if locus else '否|no'};junction={len(junctions)}|Stage 04 done")
    if locus:
        logger.info(
            f"完整locus: LB={locus.lead}bp insert={locus.insert_len}bp "
            f"RB={locus.trail}bp 总长{len(locus.seq)}bp|Complete locus built")
    paths = write_junction_outputs(junctions, contig_seqs, locus, ldir, sample)
    return {"locus": locus, "junctions": junctions, "contig_seqs": contig_seqs,
            "paths": paths, "contigs": contigs}


def run_verify_stage(cfg, logger, runner, sample: str, locus, r1: Path,
                     r2: Path, record_depth: dict, sdir: Path) -> dict:
    """05:WGS比回locus覆盖验证|05: verify coverage against rebuilt locus"""
    vdir = sdir / "05_verify"
    summary_tsv = vdir / f"{sample}.verification_summary.tsv"
    cov_tsv = vdir / f"{sample}.coverage.tsv"
    if locus is None:
        grade = "FAIL"
        with open(summary_tsv, "w") as fh:
            fh.write("sample\tgrade\tlb_junction_reads\trb_junction_reads\n"
                     f"{sample}\t{grade}\t0\t0\n")
        logger.warning("未重构出完整locus,验证分级FAIL|No locus, grade FAIL")
        return {"grade": grade, "segments": [], "lb_junction_reads": 0,
                "rb_junction_reads": 0, "summary": summary_tsv}

    locus_fa = sdir / "04_locus" / f"{sample}.complete_locus.fasta"
    vs_bam = vdir / f"{sample}.vs_locus.sorted.bam"
    if not runner.run([cfg.bwa_path, "index", str(locus_fa)],
                      description="索引locus|Index locus"):
        raise RuntimeError("locus索引失败|Index locus failed")
    ok = runner.run_pipeline(
        [[cfg.bwa_path, "mem", "-t", str(cfg.threads), str(locus_fa),
          str(r1), str(r2)],
         [cfg.samtools_path, "sort", "-@", str(cfg.threads), "-o",
          str(vs_bam), "-"]],
        description="WGS比回locus|Map WGS to locus")
    if not ok:
        raise RuntimeError("locus比对失败|Map to locus failed")
    runner.run([cfg.samtools_path, "index", str(vs_bam)],
               description="索引vs_locus|Index vs_locus")
    depth_out = runner.run_capture(
        [cfg.samtools_path, "depth", "-a", str(vs_bam)],
        description="locus深度|Locus depth")
    depth_rows = []
    if depth_out:
        for line in depth_out.splitlines():
            f = line.split("\t")
            if len(f) >= 3:
                depth_rows.append((int(f[1]), int(f[2])))
    segments = segment_coverage(depth_rows, locus)
    with open(cov_tsv, "w") as fh:
        fh.write("segment\tstart\tend\tlength\tmean_depth\tmin_depth\tzero_windows\n")
        for s in segments:
            fh.write(f"{s.name}\t{s.start}\t{s.end}\t{s.length}\t"
                     f"{s.mean_depth}\t{s.min_depth}\t{s.zero_windows}\n")
    # 跨界reads|Junction-spanning reads
    sam_out = runner.run_capture(
        [cfg.samtools_path, "view", str(vs_bam)],
        description="读vs_locus比对|Read vs_locus alignment")
    sam_lines = sam_out.splitlines() if sam_out else []
    lb_boundary = locus.lead
    rb_boundary = locus.lead + locus.insert_len
    lb_reads = count_junction_reads(sam_lines, lb_boundary)
    rb_reads = count_junction_reads(sam_lines, rb_boundary)
    flanks_plant = _check_flanks_plant(cfg, logger, runner, locus, sdir, sample)
    grade = classify(segments, lb_reads, rb_reads, cfg.junction_flank, locus,
                     flanks_plant=flanks_plant)
    with open(summary_tsv, "w") as fh:
        fh.write("sample\tgrade\tlb_junction_reads\trb_junction_reads\t"
                 "lb_depth\tinsert_depth\trb_depth\tlb_plant\trb_plant\n")
        d = {s.name: s.mean_depth for s in segments}
        lp, rp = flanks_plant if flanks_plant else ("NA", "NA")
        fh.write(f"{sample}\t{grade}\t{lb_reads}\t{rb_reads}\t"
                 f"{d.get('LB', 0)}\t{d.get('insert', 0)}\t{d.get('RB', 0)}\t"
                 f"{lp}\t{rp}\n")
    logger.info(
        f"05完成:分级={grade} LB跨界={lb_reads} RB跨界={rb_reads} "
        f"侧翼植物性={flanks_plant}|Stage 05 done")
    return {"grade": grade, "segments": segments, "lb_junction_reads": lb_reads,
            "rb_junction_reads": rb_reads, "summary": summary_tsv,
            "flanks_plant": flanks_plant}


# 侧翼匹配构建≥此bp判非植物|Flank matching construct >= this bp is non-plant
FLANK_CONSTRUCT_MATCH_BP = 50


def _check_flanks_plant(cfg, logger, runner, locus, sdir: Path, sample: str):
    """locus两端侧翼比回构建,匹配≥50bp判载体骨架来源(当年LB误判坑)|
    Align locus flanks back to construct; >=50bp match means backbone origin"""
    if locus is None:
        return None
    flank_fa = sdir / "05_verify" / f"{sample}.flank_check.fasta"
    flank_fa.parent.mkdir(parents=True, exist_ok=True)
    with open(flank_fa, "w") as fh:
        fh.write(f">lb_flank\n{locus.seq[:locus.lead]}\n")
        fh.write(f">rb_flank\n{locus.seq[len(locus.seq) - locus.trail:]}\n")
    out = runner.run_pipeline_capture(
        [[cfg.bwa_path, "mem", "-t", str(cfg.threads), str(cfg.insert_fasta),
          str(flank_fa)],
         [cfg.samtools_path, "view", "-"]],
        description="侧翼植物性校验|Flank plant check")
    if out is None:
        return None   # 校验失败不阻断分级|Check failure does not block grading
    matched = max_matched_by_query(out.splitlines())
    lb_plant = matched.get("lb_flank", 0) < FLANK_CONSTRUCT_MATCH_BP
    rb_plant = matched.get("rb_flank", 0) < FLANK_CONSTRUCT_MATCH_BP
    if not (lb_plant and rb_plant):
        logger.warning(
            f"侧翼匹配构建序列(载体骨架来源): lb_flank匹配{matched.get('lb_flank', 0)}bp "
            f"rb_flank匹配{matched.get('rb_flank', 0)}bp|"
            f"Flank(s) match construct (backbone origin)")
    return (lb_plant, rb_plant)


# ---------- 单样本/全流程|Per-sample & full pipeline ----------

def run_sample(cfg, logger, sample: str, r1: Path, r2: Path) -> dict:
    """单样本全流程(by-step输出,文件名带样本前缀)|
    Per-sample full pipeline (by-step layout, sample-prefixed files)"""
    # by-step:所有样本共享步骤目录,文件名{sample}前缀区分|
    # by-step: samples share step dirs; files distinguished by {sample} prefix
    sdir = cfg.output_path
    log_dir = sdir / "99_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    # 预建各阶段目录(samtools sort等不建父目录)|Pre-create stage dirs
    # (samtools sort cannot create parents itself)
    for stage in ("01_mapping", "02_junction_reads", "03_walking",
                  "04_locus", "05_verify"):
        (sdir / stage).mkdir(parents=True, exist_ok=True)
    runner = CommandRunner(logger, working_dir=cfg.output_path)

    # 01 比对|Mapping
    m = {}
    run_stage_with_checkpoint(
        "01_mapping",
        [sdir / "01_mapping" / f"{sample}.vs_insert.sorted.bam",
         sdir / "01_mapping" / f"{sample}.vs_insert.flagstat.txt",
         sdir / "01_mapping" / f"{sample}.insert_coverage.tsv"],
        cfg.force,
        lambda: m.update(run_mapping_stage(cfg, logger, runner, sample, r1, r2, sdir)),
        logger)
    if not m:
        cov = _read_coverage_tsv(
            sdir / "01_mapping" / f"{sample}.insert_coverage.tsv")
        bam = sdir / "01_mapping" / f"{sample}.vs_insert.sorted.bam"
        anchor_label = "insert区|insert" if "insert区|insert" in cov else None
        regions = None
        # 旧tsv缺region行(01中途崩过)→重算region统计防backbone误报|
        # Legacy tsv lacks region rows (crashed mid-stage) -> recompute
        if anchor_label is None and cfg.tdna_fasta:
            regions = _compute_insert_regions(cfg, runner, sdir)
            total_lens = read_insert_record_lens(Path(cfg.insert_fasta), sdir,
                                                 cfg.samtools_path)
            cov = _region_stats(cfg, runner, bam, sdir, regions, total_lens)
            anchor_label = "insert区|insert"
        m = {"bam": bam, "record_depth": cov, "anchor_label": anchor_label,
             "insert_regions": regions}

    # 02 junction钓取|Junction fishing
    j = {}
    run_stage_with_checkpoint(
        "02_junction_reads",
        [sdir / "02_junction_reads" / f"{sample}.softclip.fastq",
         sdir / "02_junction_reads" / f"{sample}.flank_candidates_R1.fastq",
         sdir / "02_junction_reads" / f"{sample}.flank_candidates_R2.fastq"],
        cfg.force,
        lambda: j.update(run_junction_stage(
            cfg, logger, runner, sample, m["bam"], r1, r2, sdir)),
        logger)
    if not j:
        j = {"softclip_fastq": sdir / "02_junction_reads" / f"{sample}.softclip.fastq",
             "pe1": sdir / "02_junction_reads" / f"{sample}.flank_candidates_R1.fastq",
             "pe2": sdir / "02_junction_reads" / f"{sample}.flank_candidates_R2.fastq"}

    # 03 迭代步移|Walking
    # 中间文件整体进样本子目录:master/recruited/bait等无样本前缀,
    # 混在共享目录会多样本串扰|Sample subdir: unprefixed walking
    # intermediates must not crosstalk between samples
    walk_dir = sdir / "03_walking" / "rounds" / sample
    if cfg.force and walk_dir.exists():
        # force须清掉本样本步移残留,否则旧summary导致"续跑"而非重跑|
        # force must clear this sample's walking state or the old summary
        # turns the rerun into a resume
        shutil.rmtree(walk_dir)
    walker = WalkingRunner(cfg, logger, runner)
    run_stage_with_checkpoint(
        "03_walking",
        [walk_dir / f"{sample}.walk_done.flag"],
        cfg.force,
        lambda: walker.run(sample, r1, r2, j["softclip_fastq"],
                           j["pe1"], j["pe2"], walk_dir),
        logger)

    # 04 完整locus|Complete locus
    insert_regions = m.get("insert_regions")
    if insert_regions is None and cfg.tdna_fasta:
        insert_regions = _compute_insert_regions(cfg, runner, sdir)
    l = {}
    run_stage_with_checkpoint(
        "04_locus",
        [sdir / "04_locus" / f"{sample}.junction_report.tsv",
         sdir / "04_locus" / f"{sample}.contigs.fasta"],
        cfg.force,
        lambda: l.update(run_locus_stage(
            cfg, logger, runner, sample, m["bam"], r1, r2,
            j["softclip_fastq"], walk_dir, sdir,
            insert_regions=insert_regions)),
        logger)
    if not l:
        # 续跑场景:重读产物|Resume: reload artifacts
        contigs = sdir / "04_locus" / f"{sample}.contigs.fasta"
        contig_seqs = read_fasta_dict(contigs) if contigs.exists() else {}
        contig_bam = sdir / "04_locus" / f"{sample}.contigs_vs_insert.sorted.bam"
        sam_out = runner.run_capture(
            [cfg.samtools_path, "view", str(contig_bam)],
            description="读contigs比对(续跑)|Read contig alignment (resume)")
        sam_lines = sam_out.splitlines() if sam_out else []
        if cfg.tdna_fasta:
            tdna_bam = sdir / "04_locus" / f"{sample}.contigs_vs_tdna.sorted.bam"
            tdna_sam_out = runner.run_capture(
                [cfg.samtools_path, "view", str(tdna_bam)],
                description="读contigs_vs_tdna(续跑)|Read vs tdna (resume)")
            locus_sam = tdna_sam_out.splitlines() if tdna_sam_out else []
            anchor_lens = read_insert_record_lens(Path(cfg.tdna_fasta), sdir,
                                                  cfg.samtools_path)
        else:
            locus_sam = sam_lines
            anchor_lens = read_insert_record_lens(Path(cfg.insert_fasta), sdir,
                                                  cfg.samtools_path)
        locus = find_complete_locus(locus_sam, contig_seqs, anchor_lens)
        junctions = annotate_junctions(sam_lines, contig_seqs,
                                       cfg.junction_flank, cfg.min_unmapped,
                                       insert_regions=insert_regions)
        l = {"locus": locus, "junctions": junctions, "contig_seqs": contig_seqs,
             "contigs": contigs}

    # 05 验证|Verify
    v = {}
    run_stage_with_checkpoint(
        "05_verify",
        [sdir / "05_verify" / f"{sample}.verification_summary.tsv"],
        cfg.force,
        lambda: v.update(run_verify_stage(
            cfg, logger, runner, sample, l["locus"], r1, r2,
            m.get("record_depth", {}), sdir)),
        logger)
    if not v:
        v = {"grade": "FAIL", "segments": [], "lb_junction_reads": 0,
             "rb_junction_reads": 0}

    # 06 报告|Report
    report_path = sdir / f"{sample}.insert2locus.report.html"
    run_stage_with_checkpoint(
        "06_report", [report_path], cfg.force,
        lambda: report.render_report(
            sample, l["locus"], v["segments"], v["grade"],
            v["lb_junction_reads"], v["rb_junction_reads"],
            l["junctions"], m.get("record_depth", {}), [], sdir,
            getattr(run_sample, "_versions", {}),
            junction_seqs=l["contig_seqs"],
            anchor_label=m.get("anchor_label")),
        logger)
    return {"sample": sample, "grade": v["grade"], "report": str(report_path)}


def _read_coverage_tsv(tsv: Path) -> Dict[str, dict]:
    """续跑时读回覆盖统计(跳过#注释区段行)|Reload coverage stats on resume"""
    out = {}
    if not tsv.exists():
        return out
    for line in tsv.read_text().splitlines():
        if line.startswith("#") or line.startswith("record\t"):
            continue
        f = line.split("\t")
        if len(f) >= 4:
            try:
                out[f[0]] = {"length": int(f[1]), "mean": float(f[2]),
                             "covered_frac": float(f[3])}
            except ValueError:
                continue
    return out


def main():
    """主函数|Main entry"""
    args = parse_arguments()
    try:
        cfg = Insert2locusConfig(
            input_path=args.input, insert_fasta=args.insert_fasta,
            output_dir=args.output_dir, threads=args.threads,
            sort_mem=args.sort_mem, read1_suffix=args.read1_suffix,
            max_rounds=args.max_rounds, min_softclip=args.min_softclip,
            min_unmapped=args.min_unmapped, min_growth=args.min_growth,
            mapq_min=args.mapq_min, repeat_cap=args.repeat_cap,
            junction_flank=args.junction_flank, tdna_fasta=args.tdna_fasta,
            target_flank=args.target_flank, force=args.force,
            log_level=args.log_level)
        cfg.validate()
    except ValueError as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)

    samples = cfg.discover_samples()
    info_dir = cfg.output_path / "00_pipeline_info"
    logs_dir = cfg.output_path / "99_logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    logger = Insert2locusLogger(
        logs_dir, f"insert2locus_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log",
        cfg.log_level).get_logger()
    logger.info(
        f"发现{len(samples)}个样本: {', '.join(s for s, _, _ in samples)}|"
        f"Discovered {len(samples)} sample(s)")

    # 版本信息|Version info
    runner = CommandRunner(logger, working_dir=cfg.output_path)
    versions = probe_tool_versions(cfg, runner)
    run_sample._versions = versions
    _write_pipeline_info(cfg, info_dir, versions, args)

    results, failed = [], []
    for sample, r1, r2 in samples:
        logger.info(f"===== 处理样本|Processing sample: {sample} =====")
        try:
            results.append(run_sample(cfg, logger, sample, r1, r2))
        except Exception as e:   # 单样本失败不中断其余|One failure does not stop others
            logger.error(f"样本失败|Sample failed: {sample}: {e}")
            failed.append({"sample": sample, "grade": "FAIL", "error": str(e)})

    report.render_index(results + failed, cfg.output_path / "index.html")
    # 清理临时目录|Clean tmp
    tmp_dir = cfg.output_path / "tmp"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir, ignore_errors=True)
        logger.info("已清理tmp目录|tmp cleaned")
    logger.info(
        f"完成:成功{len(results)}个,失败{len(failed)}个;导航页index.html|"
        f"Done: {len(results)} ok, {len(failed)} failed")
    sys.exit(1 if failed else 0)


def _write_pipeline_info(cfg, info_dir: Path, versions: dict, args):
    """写software_versions.yml与参数|Write software_versions.yml & params"""
    info_dir.mkdir(parents=True, exist_ok=True)
    try:
        import yaml
        info = {
            "pipeline": {"name": "biopytools insert2locus",
                         "version": __import__(
                             "biopytools.insert2locus", fromlist=["__version__"]
                         ).__version__},
            "tools": versions,
            "parameters": {
                k: getattr(args, k) for k in vars(args)
                if not k.startswith("_")},
            "execution": {"start_time": datetime.now().strftime(
                "%Y-%m-%d %H:%M:%S")},
        }
        with open(info_dir / "software_versions.yml", "w") as fh:
            yaml.dump(info, fh, default_flow_style=False, allow_unicode=True)
    except ImportError:
        info_dir.joinpath("software_versions.yml").write_text(
            "yaml模块不可用,仅记录工具版本|yaml unavailable\n" +
            "\n".join(f"{k}: {v}" for k, v in versions.items()))


if __name__ == "__main__":
    main()
