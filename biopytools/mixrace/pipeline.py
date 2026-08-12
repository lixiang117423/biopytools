"""mixrace 流程步骤(命令构造与执行)|mixrace pipeline steps (command construction & execution).

step01 run_index/run_qc, step02 run_align, step03 run_call, step04 run_filter,
step07 run_kmer。每步带断点续传;外部命令经 CommandRunner(conda 包装/直调)。
|step01-04 + 07 command runners with checkpoint resume; external commands via
CommandRunner (conda-wrapped or direct).
"""
import os
import shutil
from pathlib import Path
from typing import Optional

from .utils import get_conda_env


def run_index(config, runner, ckpt) -> Path:
    """step01a: bwa-mem2 索引参考基因组(全局,一次)|bwa-mem2 index reference (global, once)."""
    idx_dir = Path(config.output_dir) / "00_pipeline_info" / "index"
    idx_dir.mkdir(parents=True, exist_ok=True)
    fa_copy = idx_dir / os.path.basename(config.genome)
    if config.enable_checkpoint and ckpt.exists("genome_index"):
        runner.logger.info("跳过已完成步骤|Skipping completed step: genome_index")
        return idx_dir
    if not fa_copy.exists():
        shutil.copy(config.genome, str(fa_copy))   # 索引产物与 fa 同目录|index files beside fa
    runner.logger.info("开始步骤|Starting step: bwa-mem2 index")
    runner.run_conda(config.bwa_mem2_path, ["index", str(fa_copy)],
                     "bwa-mem2索引|bwa-mem2 indexing")
    if config.enable_checkpoint:
        ckpt.create("genome_index")
    return idx_dir


def run_qc(config, runner, ckpt) -> Path:
    """step01b: fastp QC 去接头(全局,处理整个 fastq 目录)|fastp QC (global, whole dir)."""
    clean_dir = Path(config.output_dir) / "01_qc"
    if config.enable_checkpoint and ckpt.exists("qc"):
        runner.logger.info("跳过已完成步骤|Skipping completed step: qc")
        return clean_dir
    clean_dir.mkdir(parents=True, exist_ok=True)
    runner.logger.info("开始步骤|Starting step: fastp QC")
    # 调兄弟命令 biopytools fastp(它内部自管 env 与配对)|sibling cmd (manages env+pairing)
    runner.run(
        f"biopytools fastp -i {config.fastq_dir} -o {clean_dir} -t {config.threads}",
        "fastp质控|fastp QC")
    if config.enable_checkpoint:
        ckpt.create("qc")
    return clean_dir


def run_align(config, runner, ckpt, sample: str, r1: str, r2: str, index_dir: str) -> Path:
    """step02: 比对 + 正确去重 + flagstat/stats|align + correct markdup + QC stats.

    关键(landmine):samtools markdup 必须 name-sort → fixmate -m → coord-sort → markdup。
    仓库 bwa/alignment.py 的实现缺 fixmate,此处正确实现。|Critical: samtools markdup
    requires name-sort → fixmate -m → coord-sort → markdup. Implemented correctly here.
    """
    out = Path(config.output_dir) / sample / "02_alignment"
    out.mkdir(parents=True, exist_ok=True)
    fa = str(Path(index_dir) / os.path.basename(config.genome))
    aln = out / f"{sample}.bam"
    name_bam = out / f"{sample}.name.bam"
    fm_bam = out / f"{sample}.fixmate.bam"
    coord_bam = out / f"{sample}.sorted.bam"
    md_bam = out / f"{sample}.sorted.markdup.bam"
    if config.enable_checkpoint and ckpt.exists(f"align_{sample}"):
        runner.logger.info(f"跳过已完成步骤|Skipping completed step: align_{sample}")
        return md_bam
    runner.logger.info(f"开始步骤|Starting step: align {sample}")
    t = config.threads
    st = config.samtools_path
    # bwa-mem2 与 samtools 同在 cphasing env → 单 conda-run bash -c 管道(禁 conda run|conda run)
    # |bwa-mem2 + samtools both in cphasing env → single conda-run bash -c pipe
    env = get_conda_env(config.bwa_mem2_path)
    # 1. bwa-mem2 mem | samtools view → 未排序 bam|unsorted bam
    runner.run(
        f'conda run -n {env} --no-capture-output bash -c '
        f'"bwa-mem2 mem -t {t} {fa} {r1} {r2} | samtools view -bS -o {aln} -"',
        f"bwa-mem2比对|bwa-mem2 mem {sample}")
    # 2. name sort(markdup 前置)|name sort (prerequisite)
    runner.run_conda(st, ["sort", "-n", "-o", str(name_bam), str(aln)],
                     f"name排序|name sort {sample}")
    # 3. fixmate -m(必须,加 mate 注释)|fixmate -m (required)
    runner.run_conda(st, ["fixmate", "-m", str(name_bam), str(fm_bam)],
                     f"fixmate|fixmate {sample}")
    # 4. coord sort|coordinate sort
    runner.run_conda(st, ["sort", "-o", str(coord_bam), str(fm_bam)],
                     f"坐标排序|coord sort {sample}")
    # 5. markdup
    runner.run_conda(st, ["markdup", "-@", str(t), str(coord_bam), str(md_bam)],
                     f"去重|markdup {sample}")
    # 6. index
    runner.run_conda(st, ["index", str(md_bam)], f"索引|index {sample}")
    # 7. flagstat → 文件|to file
    ok, flagstat_txt, _ = runner.run_conda(st, ["flagstat", str(md_bam)],
                                           f"flagstat {sample}")
    if ok and flagstat_txt:
        (out / "flagstat.txt").write_text(flagstat_txt)
    # 8. samtools stats → 文件(体量小可捕获;含 average depth)|stats to file (small, capturable)
    ok2, stats_txt, _ = runner.run_conda(st, ["stats", str(md_bam)],
                                         f"samtools stats {sample}")
    if ok2 and stats_txt:
        (out / "stats.txt").write_text(stats_txt)
    if config.enable_checkpoint:
        ckpt.create(f"align_{sample}")
    return md_bam


def parse_mean_depth(stats_text: str) -> Optional[float]:
    """从 samtools stats 输出解析平均深度|parse mean depth from samtools stats output."""
    if not stats_text:
        return None
    for line in stats_text.splitlines():
        if line.startswith("SN") and "average depth" in line:
            try:
                return float(line.split("\t")[-1])
            except (ValueError, IndexError):
                return None
    return None


def run_call(config, runner, ckpt, sample: str, bam: str, fa: str) -> Path:
    """step03: bcftools mpileup + call -mv(保留 multiallelic,带 FORMAT/AD)|variant calling.

    保留 multiallelic 位点(-m 多等位模型,不 -m2 -M2);mpileup 带 FORMAT/AD,FORMAT/DP
    供 step05 算 VAF。两 bcftools 同 sv_calling env → 单 conda-run bash -c 管道。
    |Keeps multiallelic sites (-m model, no -m2 -M2); FORMAT/AD/DP retained for VAF.
    Both bcftools in sv_calling env → single conda-run bash -c pipe.
    """
    out = Path(config.output_dir) / sample / "03_variants"
    out.mkdir(parents=True, exist_ok=True)
    raw = out / f"{sample}.raw.vcf.gz"
    if config.enable_checkpoint and ckpt.exists(f"call_{sample}"):
        runner.logger.info(f"跳过已完成步骤|Skipping completed step: call_{sample}")
        return raw
    runner.logger.info(f"开始步骤|Starting step: variant calling {sample}")
    env = get_conda_env(config.bcftools_path)
    runner.run(
        f'conda run -n {env} --no-capture-output bash -c '
        f'"bcftools mpileup -a FORMAT/AD,FORMAT/DP -f {fa} {bam} | '
        f'bcftools call -mv -Oz -o {raw}"',
        f"变异检测|variant calling {sample}")
    runner.run_conda(config.bcftools_path, ["index", "-t", str(raw)], "index vcf")
    if config.enable_checkpoint:
        ckpt.create(f"call_{sample}")
    return raw


def _filter_alt_reads(query_text: str, min_alt_reads: int):
    """从 bcftools query(CHROM/POS/AD) 筛选 max(ALT AD)>=阈值的位点|keep sites whose
    strongest ALT allele has >= min_alt_reads support. AD = ref,alt1,alt2,..."""
    keep = []
    for line in query_text.splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 3:
            continue
        chrom, pos, ad_str = f[0], f[1], f[2]
        try:
            ads = [int(x) for x in ad_str.split(",") if x and x != "."]
        except ValueError:
            continue
        alt_ads = ads[1:]                       # 首个是 ref|first is ref
        if alt_ads and max(alt_ads) >= min_alt_reads:
            keep.append(f"{chrom}\t{pos}")
    return keep


def _filter_stats_lines(counts: dict):
    """格式化过滤级联统计为 TSV 行|format filter cascade stats as TSV lines."""
    lines = ["stage\tcount"]
    for stage, n in counts.items():
        lines.append(f"{stage}\t{n}")
    return lines


def _count_records(runner, bcft: str, vcf: str) -> int:
    """统计 VCF 记录数(bcftools view -H | wc -l,流式)|count VCF records (streaming wc -l)."""
    env = get_conda_env(bcft)
    ok, out, _ = runner.run(
        f"conda run -n {env} --no-capture-output bash -c 'bcftools view -H {vcf} | wc -l'",
        "计数|count records")
    if ok and out.strip().isdigit():
        return int(out.strip())
    return 0


def run_filter(config, runner, ckpt, sample: str, raw_vcf: str) -> Path:
    """step04: VCF 过滤(QUAL/DP/SNP + 去repeat + ALT reads),保留 multiallelic|VCF filter.

    不做 -m2 -M2(保留多等位);不按 ALT reads 拆分等位,只整位点保留/丢弃。
    |No -m2 -M2 (keeps multiallelic); sites kept/dropped wholesale by ALT support.
    """
    out = Path(config.output_dir) / sample / "04_filtered"
    out.mkdir(parents=True, exist_ok=True)
    filt = out / f"{sample}.filtered.vcf.gz"
    if config.enable_checkpoint and ckpt.exists(f"filter_{sample}"):
        runner.logger.info(f"跳过已完成步骤|Skipping completed step: filter_{sample}")
        return filt
    runner.logger.info(f"开始步骤|Starting step: filter {sample}")
    bcft = config.bcftools_path
    expr = f"QUAL>={config.min_qual} && INFO/DP>={config.min_dp}"
    qd = out / "qual_dp_snp.vcf.gz"
    # QUAL + DP + 仅 SNP(保留 multiallelic SNP)|QUAL+DP+SNPs only (keeps multiallelic SNPs)
    runner.run_conda(bcft, ["view", "-v", "snps", "-i", expr, "-Oz", "-o", str(qd), str(raw_vcf)],
                     "QUAL/DP/SNP过滤|QUAL/DP/SNP filter")
    runner.run_conda(bcft, ["index", "-t", str(qd)], "index")
    cur = qd
    if config.repeat_bed:
        nr = out / "norepeat.vcf.gz"
        env = get_conda_env(bcft)
        # bcftools view | bedtools intersect -v(排除 repeat);同 sv_calling env 单 conda-run bash -c
        runner.run(
            f'conda run -n {env} --no-capture-output bash -c '
            f'"bcftools view {qd} | bedtools intersect -v -a - -b {config.repeat_bed} | '
            f'bcftools sort -Oz -o {nr}"',
            "去repeat区域|exclude repeats")
        runner.run_conda(bcft, ["index", "-t", str(nr)], "index")
        cur = nr
    # ALT reads 过滤:query AD → keep-list → view -T|ALT-reads filter via keep-list
    ok, qtxt, _ = runner.run_conda(bcft, ["query", "-f", "%CHROM\t%POS\t%AD\n", str(cur)],
                                  "query AD")
    keep = _filter_alt_reads(qtxt, config.min_alt_reads)
    keep_file = out / "keep_regions.tsv"
    keep_file.write_text("\n".join(keep) + ("\n" if keep else ""))
    runner.run_conda(bcft, ["view", "-T", str(keep_file), "-Oz", "-o", str(filt), str(cur)],
                     "ALT reads过滤|ALT-reads filter")
    runner.run_conda(bcft, ["index", "-t", str(filt)], "index")
    # 级联统计|cascade stats
    counts = {"raw": _count_records(runner, bcft, raw_vcf),
              "after_QUAL_DP_SNP": _count_records(runner, bcft, qd)}
    if config.repeat_bed:
        counts["after_repeat"] = _count_records(runner, bcft, cur)
    counts["after_ALT_reads"] = _count_records(runner, bcft, filt)
    (out / "filter_stats.tsv").write_text("\n".join(_filter_stats_lines(counts)) + "\n")
    if config.enable_checkpoint:
        ckpt.create(f"filter_{sample}")
    return filt


def run_kmer(config, runner, ckpt, clean_dir: str) -> Path:
    """step07: k-mer 谱分析(复用 biopytools smudgescope)|k-mer spectrum via smudgescope.

    smudgescope 一次处理整个 clean 目录(多样本 by-sample 输出)。读回路径:
    <kmer_root>/<sample>/02_genomescope/{model.txt,summary.txt,linear_plot.png};
    <kmer_root>/<sample>/03_smudgeplot/<sample>_smudgeplot.png。
    |smudgescope processes the whole clean dir (by-sample output). Read back from
    <root>/<sample>/02_genomescope/... and 03_smudgeplot/<sample>_smudgeplot.png.
    """
    kmer_root = Path(config.output_dir) / "07_kmer"
    if config.enable_checkpoint and ckpt.exists("kmer"):
        runner.logger.info("跳过已完成步骤|Skipping completed step: kmer")
        return kmer_root
    runner.logger.info("开始步骤|Starting step: smudgescope k-mer谱")
    runner.run(
        f"biopytools smudgescope -i {clean_dir} -o {kmer_root} "
        f"-l {config.read_length} -k {config.kmer_size} -t {config.threads}",
        "k-mer谱分析|k-mer spectrum (smudgescope)")
    if config.enable_checkpoint:
        ckpt.create("kmer")
    return kmer_root


def parse_genomescope_model(text: str) -> dict:
    """解析 genomescope2 model.txt(kcov + 杂合度 r)|parse model.txt (kcov + het r)."""
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
            elif k == "r":
                out["heterozygosity"] = v
    return out
