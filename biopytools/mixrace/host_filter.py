"""mixrace 寄主剔除|mixrace host-read depletion.

clean fastq → bwa-mem2 比对寄主 → 置信寄主(MAPQ≥Q)read 名单 → 整对剔除 → nohost fastq;
病原比对后合并统计(寄主率/mapping率/覆盖度/污染reads,即两步都未比对上的reads)。
|clean fastq -> align to host -> confident host names (MAPQ>=Q) -> drop whole pairs
-> nohost fastq; merge pathogen stats (host / map / breadth / contamination reads).
"""
import gzip
import os
import shlex
import shutil
from itertools import zip_longest
from pathlib import Path
from typing import Optional

from .pipeline import _done
from .utils import format_number, get_conda_env


def _count_bam(runner, samtools_path: str, bam, excl_flags: str,
               min_mapq: Optional[int] = None, threads: int = 1) -> Optional[int]:
    """samtools view -c 计数(失败 None)|count records (None on failure)."""
    args = ["view", "-c", "-@", str(threads), "-F", excl_flags]
    if min_mapq:
        args += ["-q", str(min_mapq)]
    args.append(str(bam))
    ok, out, _ = runner.run_conda(samtools_path, args,
                                  f"计数|count {os.path.basename(str(bam))}")
    if ok and out.strip().isdigit():
        return int(out.strip())
    return None

# samtools 过滤旗标|samtools exclude flags
_EXCL_MAPPED = "0x904"    # 排除 unmapped(0x4)+secondary(0x100)+supplementary(0x800)
NOHOST_R1_SUFFIX = "_1.nohost.fq.gz"
NOHOST_R2_SUFFIX = "_2.nohost.fq.gz"

# 统计表字段顺序(全量表)|full stats table field order
_STAT_FIELDS = ["sample", "total_reads", "host_mapped_reads", "host_rate", "nonhost_reads",
                "pathogen_mapped_reads", "pathogen_map_rate", "mean_depth", "breadth_1x",
                "unassigned_reads", "unassigned_rate"]


def _norm_read_name(header: str) -> str:
    """fastq 头 → 裸 read 名(去 @、空格注释、尾部 /1 或 /2;bwa 截掉 /1 /2 需对齐)。
    |fastq header -> bare name (strip @, comments, trailing /1 or /2; bwa trims /1 /2)."""
    name = header.lstrip("@").split()[0] if header.strip() else ""
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


def _open_reader(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def _open_writer(path: str):
    # gzip 1级压缩:中间产物喂下游GTX,速度优先(9级慢3-5x)|level-1 gzip (intermediate, speed first)
    return (gzip.open(path, "wt", compresslevel=1) if str(path).endswith(".gz")
            else open(path, "w"))


def load_host_names(names_file: str) -> set:
    """读置信寄主 read 名单(去重去空白)|load confident host read names (dedup, strip)."""
    names = set()
    p = Path(names_file)
    if not p.exists():
        return names
    with open(p, encoding="utf-8") as fh:
        for line in fh:
            n = line.strip()
            if n:
                names.add(n)
    return names


def _fq_records(fh):
    """逐条产出 4 行 fastq 记录|yield 4-line fastq records."""
    while True:
        lines = [fh.readline() for _ in range(4)]
        if not lines[0]:
            return
        if any(l == "" for l in lines[1:]):
            raise ValueError("fastq 记录不足4行(截断)|truncated fastq record (<4 lines)")
        yield [l.rstrip("\n") for l in lines]


def filter_fastq_pairs(r1_in: str, r2_in: str, r1_out: str, r2_out: str,
                       host_names: set) -> dict:
    """整对剔除:任一 mate 命中寄主名单 → 该对丢弃;输出 nohost fastq 并计数。
    |drop whole pair if either mate is a confident host read; write nohost fastq + counts."""
    total_pairs = kept_pairs = 0
    with _open_reader(r1_in) as f1, _open_reader(r2_in) as f2, \
            _open_writer(r1_out) as o1, _open_writer(r2_out) as o2:
        for rec1, rec2 in zip_longest(_fq_records(f1), _fq_records(f2)):
            if rec1 is None or rec2 is None:
                raise ValueError("R1/R2 条数不等(文件截断)|R1/R2 record counts differ "
                                 "(truncated file)")
            total_pairs += 1
            if _norm_read_name(rec1[0]) in host_names or _norm_read_name(rec2[0]) in host_names:
                continue
            kept_pairs += 1
            o1.write("\n".join(rec1) + "\n")
            o2.write("\n".join(rec2) + "\n")
    return {"total_pairs": total_pairs, "total_reads": total_pairs * 2,
            "kept_pairs": kept_pairs, "kept_reads": kept_pairs * 2,
            "dropped_pairs": total_pairs - kept_pairs}


def _safe_div(a, b) -> float:
    return a / b if b else 0.0


def compute_host_stats(total_reads: int, host_mapped: Optional[int],
                       nonhost_reads: Optional[int], pathogen_mapped: int) -> dict:
    """汇总口径:寄主率/mapping率/未归属。host_mapped=None 表示未启用寄主剔除。
    |accounting: host rate / pathogen map rate / unassigned. host_mapped=None = no host stage."""
    if host_mapped is None:
        host_rate = None
        nonhost = total_reads
    else:
        host_rate = _safe_div(host_mapped, total_reads)
        nonhost = nonhost_reads if nonhost_reads is not None else total_reads - host_mapped
    unassigned = max(nonhost - pathogen_mapped, 0)
    return {
        "total_reads": total_reads, "host_mapped_reads": host_mapped,
        "host_rate": host_rate, "nonhost_reads": nonhost,
        "pathogen_mapped_reads": pathogen_mapped,
        "pathogen_map_rate": _safe_div(pathogen_mapped, nonhost),
        "unassigned_reads": unassigned,
        "unassigned_rate": _safe_div(unassigned, total_reads),
    }


def parse_samtools_coverage(text: str) -> dict:
    """聚合 samtools coverage 输出(covbases 列按表头定位)|aggregate coverage output.

    列(samtools 1.23)|columns: #rname startpos endpos numreads covbases coverage meandepth ...
    """
    cov_idx = 4
    total = 0
    contigs = 0
    for line in (text or "").splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        if line.startswith("#"):
            try:
                cov_idx = f.index("covbases")
            except ValueError:
                pass
            continue
        if len(f) > cov_idx:
            try:
                total += int(float(f[cov_idx]))
                contigs += 1
            except ValueError:
                continue
    return {"covbases_total": total, "contigs": contigs}


def build_host_stats_tsv(stats: dict) -> str:
    """全量统计表(field\tvalue;None → NA)|full stats table (None -> NA)."""
    lines = ["field\tvalue"]
    for k in _STAT_FIELDS:
        if k in stats:
            v = stats[k]
            lines.append(f"{k}\t{'NA' if v is None else v}")
    return "\n".join(lines) + "\n"


def _read_field_tsv(path) -> dict:
    """读 field\tvalue 表,数值自动转型|read field/value table, auto-convert numerics."""
    out: dict = {}
    p = Path(path)
    if not p.exists():
        return out
    for line in p.read_text(encoding="utf-8").splitlines():
        f = line.split("\t")
        if len(f) < 2 or f[0] == "field":
            continue
        v: object = f[1]
        for cast in (int, float):
            try:
                v = cast(f[1])
                break
            except ValueError:
                continue
        out[f[0]] = v
    return out


def read_mapq_stats(config, sample: str) -> dict:
    """读比对阶段计数表|read align-stage counts.

    v0.3 写在 04_het_eval/{sample}.mapq_stats.tsv;02_alignment 为 v0.2 遗留回退。
    |v0.3 path first; 02_alignment kept as v0.2 legacy fallback.
    """
    for sub in ("04_het_eval", "02_alignment"):
        d = _read_field_tsv(Path(config.output_dir) / sub / f"{sample}.mapq_stats.tsv")
        if d:
            return d
    return {}


def run_host_index(config, runner, ckpt) -> Path:
    """寄主基因组 bwa-mem2 索引(00_pipeline_info/index_host)|bwa-mem2 index for host genome."""
    idx_dir = Path(config.output_dir) / "00_pipeline_info" / "index_host"
    idx_dir.mkdir(parents=True, exist_ok=True)
    fa_copy = idx_dir / os.path.basename(config.host_genome)
    if config.enable_checkpoint and _done(ckpt, "host_index", str(fa_copy) + ".0123"):
        runner.logger.info("跳过已完成步骤|Skipping completed step: host_index")
        return idx_dir
    if not fa_copy.exists():
        shutil.copy(config.host_genome, str(fa_copy))
    runner.logger.info("开始步骤|Starting step: bwa-mem2 index (host)")
    ok, _, _ = runner.run_conda(config.bwa_mem2_path, ["index", str(fa_copy)],
                                "寄主bwa-mem2索引|host bwa-mem2 indexing")
    if config.enable_checkpoint and ok:
        ckpt.create("host_index")
    elif config.enable_checkpoint:
        runner.logger.error("host_index 失败,未建断点|host_index failed, no checkpoint")
    return idx_dir


def run_host_filter(config, runner, ckpt, sample: str, r1: str, r2: str,
                    host_index_dir: str) -> Optional[dict]:
    """逐样本寄主剔除:比对→名单→整对剔除→nohost fastq+阶段统计表。
    |per-sample host depletion: align -> names -> pair-drop -> nohost fastq + stage stats."""
    out_dir = Path(config.output_dir) / "02_host_filter"
    tmp_dir = Path(config.output_dir) / "tmp"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    nohost_r1 = out_dir / f"{sample}{NOHOST_R1_SUFFIX}"
    nohost_r2 = out_dir / f"{sample}{NOHOST_R2_SUFFIX}"
    partial_tsv = out_dir / f"{sample}.host_filter.tsv"
    if config.enable_checkpoint and _done(ckpt, f"host_filter_{sample}", nohost_r1):
        runner.logger.info(f"跳过已完成步骤|Skipping completed step: host_filter_{sample}")
        stats = _read_field_tsv(partial_tsv)
        return {**stats, "sample": sample, "nohost_r1": str(nohost_r1),
                "nohost_r2": str(nohost_r2)}

    fa = str(Path(host_index_dir) / os.path.basename(config.host_genome))
    host_bam = tmp_dir / f"{sample}.host.bam"
    names_file = tmp_dir / f"{sample}.host.names"
    runner.logger.info(f"开始步骤|Starting step: host filter {sample}")
    t = config.threads
    q = config.min_mapq
    env = get_conda_env(config.bwa_mem2_path)   # bwa-mem2 的 conda 环境
    # §13.2.3: 管道内一律完整路径(shlex.quote), 不依赖 env PATH 解析命令名
    # |full quoted paths inside the pipe (no bare command names, §13.2.3)
    bwa_path = shlex.quote(config.bwa_mem2_path)
    st_path = shlex.quote(config.samtools_path)
    # ① 寄主比对(临时 BAM,bwa 输出天然按名成对,免排序)|align to host (temp BAM)
    ok_aln, _, _ = runner.run(
        f'conda run -n {env} --no-capture-output bash -c '
        f'"{bwa_path} mem -t {t} {fa} {r1} {r2} | {st_path} view -b -@ {t} -o {host_bam} -"',
        f"寄主比对|host align {sample}")
    if not ok_aln:
        runner.logger.error(f"寄主比对失败,中止 host_filter {sample}|host align failed, aborted")
        return None
    if getattr(runner, "dry_run", False):
        # dry-run:命令已记录,不执行计数/过滤;返回零统计(断点自愈,输出缺失不会误跳过)
        # |dry-run: commands logged only; zero stats (self-healing checkpoint prevents skips)
        runner.logger.info(f"dry-run:跳过寄主计数与fastq过滤|dry-run: skip count/filter {sample}")
        return {"sample": sample, "total_reads": 0, "host_mapped_reads": 0,
                "host_rate": 0.0, "nonhost_reads": 0,
                "nohost_r1": str(nohost_r1), "nohost_r2": str(nohost_r2)}
    # ② 计数(总 primary / 置信寄主 MAPQ≥Q)|counts (total primary / confident host)
    total_primary = _count_bam(runner, config.samtools_path, str(host_bam), _EXCL_MAPPED,
                               threads=config.threads)
    host_mapped = _count_bam(runner, config.samtools_path, str(host_bam), _EXCL_MAPPED,
                             min_mapq=q, threads=config.threads)
    if total_primary is None or host_mapped is None:
        runner.logger.error(f"寄主计数失败,中止 host_filter {sample}|host count failed, aborted")
        return None
    # ③ 置信寄主 read 名单|extract confident host read names
    st_env = get_conda_env(config.samtools_path)
    ok_names, _, _ = runner.run(
        f'conda run -n {st_env} --no-capture-output bash -c '
        f"'{st_path} view -@ {t} -F {_EXCL_MAPPED} -q {q} {host_bam} | cut -f1 | sort -u --parallel={t} > {names_file}'",
        f"提取寄主read名|extract host read names {sample}")
    if not (ok_names and names_file.exists()):
        runner.logger.error(f"寄主名单提取失败,中止 host_filter {sample}|name extraction failed")
        return None
    host_names = load_host_names(str(names_file))
    # ④ 整对剔除|drop host pairs
    try:
        res = filter_fastq_pairs(r1, r2, str(nohost_r1), str(nohost_r2), host_names)
    except (OSError, ValueError) as e:
        runner.logger.error(f"fastq 过滤失败|fastq filtering failed {sample}: {e}")
        for p in (nohost_r1, nohost_r2):   # 清残留半成品|drop partial outputs
            try:
                p.unlink()
            except OSError:
                pass
        return None
    # ⑤ 清理临时文件(省空间)|clean temp
    for p in (host_bam, names_file):
        try:
            p.unlink()
        except OSError:
            pass
    total_reads = res["total_reads"]
    host_rate = _safe_div(host_mapped, total_reads)
    runner.logger.info(
        f"{sample}: 总reads|total reads: {format_number(total_reads)}, "
        f"置信寄主|confident host: {format_number(host_mapped)} ({host_rate * 100:.2f}%), "
        f"剔除非寄主对|kept non-host pairs: {res['kept_pairs']}")
    stats = {"sample": sample, "total_reads": total_reads,
             "host_mapped_reads": host_mapped, "host_rate": host_rate,
             "nonhost_reads": res["kept_reads"]}
    partial_tsv.write_text(build_host_stats_tsv(stats), encoding="utf-8")
    if config.enable_checkpoint:
        ckpt.create(f"host_filter_{sample}")
    return {**stats, "nohost_r1": str(nohost_r1), "nohost_r2": str(nohost_r2)}


def pathogen_alignment_stats(config, runner, sample: str, bam: str,
                             genome_size: int, mean_depth=None) -> dict:
    """合并寄主阶段+病原比对统计,写全量表 02_host_filter/{sample}.host_stats.tsv。

    寄主率来自 host_filter 阶段表;病原 mapped/总数来自 align 阶段 mapq_stats.tsv
    (缺失时对最终 BAM 现场计数,过滤后 BAM 无 unmapped,总数退化为 mapped 并告警)。
    |merge host-stage + pathogen-align stats into full table. Pathogen counts from
    align-stage mapq_stats.tsv; fallback = live counts on final BAM (degraded, warned).
    """
    out_dir = Path(config.output_dir) / "02_host_filter"
    out_dir.mkdir(parents=True, exist_ok=True)
    host = _read_field_tsv(out_dir / f"{sample}.host_filter.tsv")
    counts = read_mapq_stats(config, sample)
    pathogen_mapped = counts.get("mapped_q_reads")
    if host:
        total_reads = host.get("total_reads", 0)
        host_mapped = host.get("host_mapped_reads")
        nonhost = host.get("nonhost_reads")
    else:
        total_reads = counts.get("total_primary_reads")
        host_mapped = None
        nonhost = None
    if pathogen_mapped is None:
        runner.logger.warning(
            f"{sample}: mapq_stats.tsv 缺失,现场计数最终 BAM(过滤后 BAM 无 unmapped,总数退化)"
            f"|mapq_stats.tsv missing, live counts on final BAM (degraded total)")
        pathogen_mapped = _count_bam(runner, config.samtools_path, bam, _EXCL_MAPPED,
                                     threads=config.threads) or 0
        if total_reads is None:
            total_reads = pathogen_mapped
    # 覆盖广度(≥1x;samtools coverage -q 与 MAPQ 统一口径)|breadth >=1x with min-MQ
    cov_args = ["coverage"]
    if getattr(config, "min_mapq", 0) > 0:
        cov_args += ["-q", str(config.min_mapq)]
    cov_args.append(bam)
    ok_cov, cov_text, _ = runner.run_conda(config.samtools_path, cov_args,
                                           f"覆盖广度|breadth {sample}")
    cov = parse_samtools_coverage(cov_text if ok_cov else "")
    breadth = cov["covbases_total"] / genome_size * 100 if genome_size else 0.0
    stats = compute_host_stats(total_reads or 0, host_mapped, nonhost, pathogen_mapped)
    stats.update({"sample": sample, "mean_depth": mean_depth, "breadth_1x": breadth})
    (out_dir / f"{sample}.host_stats.tsv").write_text(build_host_stats_tsv(stats),
                                                      encoding="utf-8")
    runner.logger.info(
        f"{sample}: 病原mapping率|pathogen map rate: {stats['pathogen_map_rate'] * 100:.2f}%, "
        f"覆盖广度|breadth>=1x: {breadth:.2f}%, "
        f"污染reads|contamination: {format_number(stats['unassigned_reads'])} "
        f"({stats['unassigned_rate'] * 100:.2f}%)")
    return stats
