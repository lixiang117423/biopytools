"""mixrace 流程步骤(单倍体 freebayes 后端)|mixrace pipeline steps (haploid freebayes backend).

step01 run_index/run_qc(bwa-mem2 索引 + fastp QC),
step02 run_align(bwa-mem2 + 正确 markdup 五步),
step03 run_call_freebayes(freebayes -p 1 单倍体,保低频;AF 用 AO/RO),
step04 run_filter(QUAL/DP/去repeat,保留多等位),
step06 run_kmer(smudgescope)。
|Haploid pipeline: bwa-mem2+markdup, freebayes -p 1 (keeps low-freq alleles;
AF via AO/RO), filter (keeps multiallelic), smudgescope.
"""
import os
import shutil
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


def run_align(config, runner, ckpt, sample: str, r1: str, r2: str, index_dir: str) -> Path:
    """step02: bwa-mem2 比对 + 正确 markdup(name-sort→fixmate -m→sort→markdup→index)+ flagstat/stats。
    |bwa-mem2 align + correct markdup + QC. Landmine: fixmate -m required before markdup."""
    out = Path(config.output_dir) / "02_alignment"
    out.mkdir(parents=True, exist_ok=True)
    fa = str(Path(index_dir) / os.path.basename(config.genome))
    aln = out / f"{sample}.bam"
    name_bam = out / f"{sample}.name.bam"
    fm_bam = out / f"{sample}.fixmate.bam"
    coord_bam = out / f"{sample}.sorted.bam"
    md_bam = out / f"{sample}.sorted.markdup.bam"
    if config.enable_checkpoint and _done(ckpt, f"align_{sample}", md_bam):
        runner.logger.info(f"跳过已完成步骤|Skipping completed step: align_{sample}")
        return md_bam
    runner.logger.info(f"开始步骤|Starting step: align {sample}")
    t = config.threads
    st = config.samtools_path
    env = get_conda_env(config.bwa_mem2_path)   # cphasing(bwa-mem2+samtools 同 env)
    # 顺序执行,任一步失败立即短路(否则对不存在的 aln 连环报错,淹没真正的失败原因)|
    # sequential with short-circuit: later steps on a missing aln just drown the real error.
    steps = [
        (lambda: runner.run(
            # -R 读组:freebayes 从 @RG SM 取样品名;不加则所有 VCF 样品名都是 "unknown",
            # bcftools merge 报 Duplicate sample names
            # |read group: freebayes takes sample name from @RG SM; without it all
            # VCFs are named "unknown" -> merge fails with Duplicate sample names
            f'conda run -n {env} --no-capture-output bash -c '
            f'"bwa-mem2 mem -t {t} -R \'@RG\\tID:{sample}\\tSM:{sample}\' {fa} {r1} {r2} '
            f'| samtools view -bS -o {aln} -"',
            f"bwa-mem2比对|bwa-mem2 mem {sample}"), "bwa-mem2 mem"),
        (lambda: runner.run_conda(st, ["sort", "-n", "-o", str(name_bam), str(aln)],
                                  f"name排序|name sort {sample}"), "samtools sort -n"),
        (lambda: runner.run_conda(st, ["fixmate", "-m", str(name_bam), str(fm_bam)],
                                  f"fixmate|fixmate {sample}"), "samtools fixmate"),
        (lambda: runner.run_conda(st, ["sort", "-o", str(coord_bam), str(fm_bam)],
                                  f"坐标排序|coord sort {sample}"), "samtools sort"),
        (lambda: runner.run_conda(st, ["markdup", "-@", str(t), str(coord_bam), str(md_bam)],
                                  f"去重|markdup {sample}"), "samtools markdup"),
        (lambda: runner.run_conda(st, ["index", str(md_bam)], f"索引|index {sample}"),
         "samtools index"),
    ]
    bam_ok = True
    for fn, step_name in steps:
        ok, _, _ = fn()
        if not ok:
            runner.logger.error(
                f"align {sample} 步骤 [{step_name}] 失败,中断后续步骤(未建断点,重跑将重试)"
                f"|align step [{step_name}] failed for {sample}, aborting rest "
                f"(no checkpoint; rerun will retry)")
            bam_ok = False
            break
    # 清理中间 bam(省空间,by_step 目录整洁)|clean intermediate bams
    for tmp in (aln, name_bam, fm_bam, coord_bam):
        try:
            tmp.unlink()
        except OSError:
            pass
    if config.enable_checkpoint and bam_ok:
        ckpt.create(f"align_{sample}")
    elif config.enable_checkpoint:
        runner.logger.error(f"align {sample} 失败,未建断点|align failed, no checkpoint")
    return md_bam


def run_call_freebayes(config, runner, ckpt, sample: str, bam: str, fa: str) -> Path:
    """step03: freebayes 单倍体 calling(-p 1,保低频 --min-alternate-fraction)|haploid freebayes.
    AF 用 AO/RO 字段(非 AD)。freebayes|bgzip 同 freebayes env;index 走 bcftools(sv_calling)。"""
    out = Path(config.output_dir) / "03_variants"
    out.mkdir(parents=True, exist_ok=True)
    vcf = out / f"{sample}.raw.vcf.gz"
    if config.enable_checkpoint and _done(ckpt, f"call_{sample}", vcf):
        runner.logger.info(f"跳过已完成步骤|Skipping completed step: call_{sample}")
        return vcf
    runner.logger.info(f"开始步骤|Starting step: freebayes -p 1 {sample}")
    env = get_conda_env(config.freebayes_path)   # freebayes(同时有 bgzip)
    af = config.freebayes_min_alternate_fraction
    cov = config.freebayes_min_coverage
    o1, _, _ = runner.run(
        f'conda run -n {env} --no-capture-output bash -c '
        f'"freebayes -p 1 --min-alternate-fraction {af} --min-coverage {cov} '
        f'-f {fa} {bam} | bgzip -c > {vcf}"',
        f"freebayes 单倍体 calling|freebayes -p 1 {sample}")
    o2, _, _ = runner.run_conda(config.bcftools_path, ["index", "-t", str(vcf)], "index vcf")
    if config.enable_checkpoint and o1 and o2:
        ckpt.create(f"call_{sample}")
    elif config.enable_checkpoint:
        runner.logger.error(f"call {sample} 失败,未建断点|call failed, no checkpoint")
    return vcf


def _filter_stats_lines(counts: dict):
    lines = ["stage\tcount"]
    for stage, n in counts.items():
        lines.append(f"{stage}\t{n}")
    return lines


def _count_records(runner, bcft: str, vcf: str) -> int:
    env = get_conda_env(bcft)
    ok, out, _ = runner.run(
        f"conda run -n {env} --no-capture-output bash -c 'bcftools view -H {vcf} | wc -l'",
        "计数|count records")
    if ok and out.strip().isdigit():
        return int(out.strip())
    return 0


def run_filter(config, runner, ckpt, sample: str, raw_vcf: str) -> Path:
    """step04: 过滤(QUAL/DP/SNP + 去repeat),保留多等位|filter (keeps multiallelic)."""
    out = Path(config.output_dir) / "04_filtered"
    out.mkdir(parents=True, exist_ok=True)
    filt = out / f"{sample}.filtered.vcf.gz"
    if config.enable_checkpoint and _done(ckpt, f"filter_{sample}", filt):
        runner.logger.info(f"跳过已完成步骤|Skipping completed step: filter_{sample}")
        return filt
    runner.logger.info(f"开始步骤|Starting step: filter {sample}")
    bcft = config.bcftools_path
    expr = f"QUAL>={config.min_qual} && INFO/DP>={config.min_dp}"
    cur = raw_vcf
    if config.repeat_bed:
        nr = out / "norepeat.vcf.gz"
        env = get_conda_env(bcft)
        ok_r, _, _ = runner.run(
            f'conda run -n {env} --no-capture-output bash -c '
            f'"bcftools view -i \'{expr}\' {raw_vcf} | bedtools intersect -v -a - -b {config.repeat_bed} | '
            f'bcftools sort -Oz -o {nr}"', "QUAL/DP+去repeat|QUAL/DP+exclude repeats")
        if not (ok_r and nr.exists()):
            # 失败即中止:严禁回退 raw(否则 filt 未过滤却建断点,下游全基于脏数据)|
            # Abort on failure: never fall back to raw (unfiltered filt + checkpoint poisons downstream).
            runner.logger.error(
                f"去repeat管道失败,中止 filter {sample}(不写输出/不建断点,重跑将重试)"
                f"|Repeat-exclude pipeline failed, aborting filter {sample} "
                f"(no output/no checkpoint; rerun will retry)")
            try:
                filt.unlink()   # 清掉可能残留的旧 filt,防下游误用|drop stale filt
            except OSError:
                pass
            return filt
        cur = str(nr)
        runner.run_conda(bcft, ["view", "-Oz", "-o", str(filt), str(cur)], "收尾|finalize")
    else:
        runner.run_conda(bcft, ["view", "-i", expr, "-Oz", "-o", str(filt), str(raw_vcf)],
                         "QUAL/DP过滤|QUAL/DP filter")
    ok_idx, _, _ = runner.run_conda(bcft, ["index", "-t", str(filt)], "index")
    counts = {"raw": _count_records(runner, bcft, raw_vcf), "filtered": _count_records(runner, bcft, filt)}
    (out / "filter_stats.tsv").write_text("\n".join(_filter_stats_lines(counts)) + "\n")
    if config.enable_checkpoint and ok_idx and filt.exists():
        ckpt.create(f"filter_{sample}")
    elif config.enable_checkpoint:
        runner.logger.error(f"filter {sample} 失败,未建断点|filter failed, no checkpoint")
    return filt


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
    """samtools stats → alignment_qc/{sample}.stats.txt;平均深度 = bases_mapped / genome_size。
    |samtools stats -> stats.txt; mean depth = bases_mapped / genome_size."""
    qc_dir = Path(config.output_dir) / "alignment_qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    ok, stats_txt, _ = runner.run_conda(config.samtools_path, ["stats", str(bam)],
                                        f"samtools stats {sample}")
    if not (ok and stats_txt):
        return None
    (qc_dir / f"{sample}.stats.txt").write_text(stats_txt)
    bases = _parse_sn(stats_txt, "bases mapped")
    if bases and genome_size:
        return bases / genome_size
    return parse_mean_depth(stats_txt)   # 回退:某些 stats 有 average depth 行|fallback


def run_kmer(config, runner, ckpt, clean_dir: str) -> Path:
    """step06: k-mer 谱(smudgescope,读 clean fastq)|k-mer spectrum via smudgescope."""
    kmer_root = Path(config.output_dir) / "06_kmer"
    if config.enable_checkpoint and _done(ckpt, "kmer", kmer_root):
        runner.logger.info("跳过已完成步骤|Skipping completed step: kmer")
        return kmer_root
    runner.logger.info("开始步骤|Starting step: smudgescope k-mer谱")
    ok, _, _ = runner.run(
        f"biopytools smudgescope -i {clean_dir} -o {kmer_root} "
        f"-l {config.read_length} -k {config.kmer_size} -t {config.threads}",
        "k-mer谱分析|k-mer spectrum (smudgescope)")
    if config.enable_checkpoint and ok:
        ckpt.create("kmer")
    elif config.enable_checkpoint:
        runner.logger.error("kmer 失败,未建断点|kmer failed, no checkpoint")
    return kmer_root


def run_tree(config, runner, ckpt, samples: list):
    """step07b: 样品系统发育树(聚类)|phylogenetic tree for sample clustering.

    bcftools merge 各样品 filtered VCF → biopytools vcf2tree(IQ-TREE2)→ newick。
    返回 nwk 路径或 None(样品<4 / VCF 缺 / --skip-tree);HTML 用 phylocanvas.gl 交互渲染。
    |merge per-sample filtered VCFs -> vcf2tree (IQ-TREE2) -> newick.
    Returns nwk path or None (<4 samples / missing VCF / --skip-tree);
    HTML renders it interactively via phylocanvas.gl.
    """
    if config.skip_tree:
        runner.logger.info("跳过建树(--skip-tree)|skipping tree (--skip-tree)")
        return None
    tree_dir = Path(config.output_dir) / "08_tree"
    nwk = tree_dir / "vcf2tree" / "02_tree" / "merged.iqtree.nwk"
    if config.enable_checkpoint and _done(ckpt, "tree", nwk):
        runner.logger.info("跳过已完成步骤|Skipping completed step: tree")
        return nwk
    # 收集存在的 filtered VCF|collect existing filtered VCFs
    vcfs = []
    for s in samples:
        name = s["sample"] if isinstance(s, dict) else s
        v = Path(config.output_dir) / "04_filtered" / f"{name}.filtered.vcf.gz"
        if v.exists():
            vcfs.append(str(v))
        else:
            runner.logger.warning(f"建树缺 VCF,跳过该样品|tree: missing VCF, skipped: {name}")
    if len(vcfs) < 4:
        runner.logger.warning(f"可用样品 {len(vcfs)}<4,不建树(树需 >=4 样本)"
                              f"|{len(vcfs)} usable samples <4, tree skipped")
        return None
    tree_dir.mkdir(parents=True, exist_ok=True)
    merged = tree_dir / "merged.vcf.gz"
    # 0. reheader:旧 bam(无 @RG)产出的 VCF 样品名全是 "unknown",merge 会报 Duplicate
    # sample names → 每样品重命名成真实样名(幂等:已正确命名的 VCF 重写同名无副作用)
    # |reheader: VCFs from RG-less bams are all named "unknown" (merge fails on duplicates)
    # -> rename each to its real sample name (idempotent for already-named VCFs)
    bcft = config.bcftools_path
    renamed = []
    for s in samples:
        name = s["sample"] if isinstance(s, dict) else s
        src = Path(config.output_dir) / "04_filtered" / f"{name}.filtered.vcf.gz"
        if not src.exists():
            continue
        dst = tree_dir / f"{name}.renamed.vcf.gz"
        sfile = tree_dir / f"{name}.samples.txt"
        sfile.write_text(name + "\n")
        ok_r, _, _ = runner.run_conda(
            bcft, ["reheader", "-s", str(sfile), "-o", str(dst), str(src)],
            f"重命名样品|reheader {name}")
        ok_ri, _, _ = runner.run_conda(bcft, ["index", "-t", str(dst)], f"index {name}")
        if ok_r and ok_ri and dst.exists():
            renamed.append(str(dst))
        else:
            runner.logger.warning(f"reheader 失败,建树跳过该样品|reheader failed, skipped: {name}")
    if len(renamed) < 4:
        runner.logger.warning(f"可合并样品 {len(renamed)}<4,不建树|{len(renamed)} samples <4, tree skipped")
        return None
    # 1. bcftools merge(多样本合并)|merge samples
    ok_m, _, _ = runner.run_conda(
        bcft,
        ["merge", "-Oz", "-o", str(merged)] + renamed,
        f"合并 {len(renamed)} 样品 VCF|merge {len(renamed)} VCFs")
    if not ok_m:
        # 失败残壳清理,避免留下非 BGZF 空文件|clean failed stub
        if merged.exists():
            try:
                merged.unlink()
            except OSError:
                pass
        runner.logger.error("merge 失败,不建断点|merge failed, no checkpoint")
        return None
    ok_i, _, _ = runner.run_conda(bcft, ["index", "-t", str(merged)], "index merged")
    if not (ok_i and merged.exists()):
        runner.logger.error("index merged 失败,不建断点|index merged failed, no checkpoint")
        return None
    # 2. vcf2tree(IQ-TREE2,复用兄弟模块)|build tree via sibling module
    v2t_out = tree_dir / "vcf2tree"
    ok_t, _, _ = runner.run(
        f"biopytools vcf2tree -i {merged} -o {v2t_out} -t {config.threads}",
        "系统发育树(IQ-TREE2)|phylogenetic tree (IQ-TREE2)")
    if not (ok_t and nwk.exists()):
        runner.logger.error("vcf2tree 未产出 newick,不建断点|vcf2tree produced no newick, no checkpoint")
        return None
    if config.enable_checkpoint:
        ckpt.create("tree")
    return nwk


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
