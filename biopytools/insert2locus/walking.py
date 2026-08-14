"""迭代步移引擎|Iterative walking engine

把run.sh第6节逻辑Python化:侧翼诱饵→招募新reads→SPAdes重组装→循环至收敛|
Python port of run.sh section 6: flank bait -> recruit reads -> reassemble -> loop
"""

import shutil
import tempfile
import unicodedata
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

from .data_processing import extract_flank_bait, mate_suffix_names, read_fasta_dict
from .locus_builder import max_flank_extent
from .utils import build_conda_command, format_number

STOP_REASONS = {
    "no_bait": "无侧翼诱饵|No flank bait",
    "no_new_reads": "无新read,收敛|No new reads, converged",
    "repeat_cap": "新招reads超上限,疑似撞重复区,保留上一轮|Exceeded repeat cap, "
                  "suspected repeat region, keeping last round",
    "no_growth": "contig总长增量低于阈值,收敛|Growth below threshold, converged",
    "max_rounds": "达到最大轮数|Reached max rounds",
}


def decide_stop(new_reads: int, growth: int, has_bait: bool, cfg,
                current_flank: Optional[int] = None) -> Optional[str]:
    """判定是否停止及原因|Decide whether to stop and why

    侧翼未达target_flank时,小增量(<min_growth)不触发收敛,继续向目标走;
    零增长/无新reads/撞重复区仍硬停|max_rounds兜底|
    While below target_flank, small growth does not trigger convergence;
    zero growth / no new reads / repeat cap still stop hard; max_rounds caps
    """
    if not has_bait:
        return "no_bait"
    if new_reads == 0:
        return "no_new_reads"
    if new_reads > cfg.repeat_cap:
        return "repeat_cap"
    if growth < cfg.min_growth:
        if current_flank is not None and 0 < growth and current_flank < cfg.target_flank:
            return None   # 未达标仍在生长,继续|Still growing toward target
        return "no_growth"
    return None


def contig_stats(fasta: Path) -> Tuple[int, int, int]:
    """纯Python统计(total,longest,n)|Pure-Python stats (total,longest,count)"""
    total = longest = n = 0
    cur = 0
    for line in Path(fasta).read_text().splitlines():
        if line.startswith(">"):
            if n > 0:
                total += cur
                longest = max(longest, cur)
            n += 1
            cur = 0
        else:
            cur += len(line.strip())
    if n > 0:
        total += cur
        longest = max(longest, cur)
    return (total, longest, n)


def parse_summary_tsv(path: Path) -> List[dict]:
    """读回步移汇总(断点续跑用)|Read back walking summary (for resume)"""
    rows = []
    text = Path(path).read_text()
    lines = text.splitlines()
    if not lines:
        return rows
    header = lines[0].split("\t")
    for line in lines[1:]:
        vals = line.split("\t")
        rows.append({h: _num(v) for h, v in zip(header, vals)})
    return rows


def _num(v: str):
    return int(v) if v.lstrip("-").isdigit() else v


def spades_work_dir(output_tmp: Path) -> Path:
    """SPAdes工作目录:ASCII优先output_dir/tmp,路径含中文回退系统/tmp|
    SPAdes work dir: output_dir/tmp if ASCII, else system /tmp fallback
    (SPAdes对非ASCII路径不稳定|SPAdes is unstable on non-ASCII paths)"""
    if str(output_tmp).isascii():
        return Path(output_tmp)
    return Path(tempfile.mkdtemp(prefix="insert2locus_spades."))


def fish_mates(seqkit_path: str, runner, names: List[str], mate: int,
               fastq: Path, pattern_file: Path, out_fastq: Path,
               append: bool) -> bool:
    """按read名从原始fastq钓mate|Fish mates from original fastq by read name"""
    pattern_file.parent.mkdir(parents=True, exist_ok=True)
    pattern_file.write_text("\n".join(mate_suffix_names(names, mate)) + "\n")
    out = runner.run_capture(
        [seqkit_path, "grep", "-f", str(pattern_file), str(fastq)],
        description=f"seqkit钓取mate/{mate}|seqkit fish mate/{mate}")
    if out is None:
        return False
    out_fastq.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if append else "w"
    with open(out_fastq, mode) as fh:
        fh.write(out)
        if out and not out.endswith("\n"):
            fh.write("\n")
    return True


@dataclass
class WalkResult:
    """步移结果|Walking result"""
    final_contigs: Path
    rounds_done: int
    stop_reason: str          # STOP_REASONS的key|Key of STOP_REASONS
    total_bp: int
    longest_bp: int


class WalkingRunner:
    """步移编排器|Walking orchestrator"""

    def __init__(self, cfg, logger, runner):
        self.cfg = cfg
        self.logger = logger
        self.runner = runner

    # ---------- 基础步骤|Primitive steps ----------

    def run_spades(self, sc_fq: Path, pe1: Path, pe2: Path,
                   out_contigs: Path) -> bool:
        """跑SPAdes(conda run,ASCII工作目录)|Run SPAdes (conda run, ASCII work dir)"""
        out_contigs.parent.mkdir(parents=True, exist_ok=True)
        tmp_root = self.cfg.output_path / "tmp"
        tmp_root.mkdir(parents=True, exist_ok=True)
        work = spades_work_dir(tmp_root)
        work = Path(work) if work == tmp_root else work
        if work == tmp_root:
            # ASCII场景用子目录隔离|Isolate in a subdir in the ASCII case
            work = Path(tempfile.mkdtemp(prefix="spades_round.", dir=str(work)))
        try:
            w_sc = work / "softclip.fastq"
            w_p1 = work / "pe1_1.fastq"
            w_p2 = work / "pe1_2.fastq"
            shutil.copy(str(sc_fq), str(w_sc))
            shutil.copy(str(pe1), str(w_p1))
            shutil.copy(str(pe2), str(w_p2))
            args = ["-s", str(w_sc), "--pe1-1", str(w_p1), "--pe1-2", str(w_p2),
                    "-o", str(work / "out"), "--careful",
                    "-t", str(self.cfg.threads), "-m", "800"]
            cmd = build_conda_command(self.cfg.spades_path, args)
            spades_log = work / "spades.log"
            if not self.runner.run(cmd, description=f"SPAdes局部组装|SPAdes local assembly"):
                return False
            _ = spades_log  # spades输出已由runner记stderr;日志文件留存work内
            contigs = work / "out" / "contigs.fasta"
            if not contigs.exists():
                self.logger.error(
                    f"SPAdes未产出contigs|SPAdes produced no contigs: {contigs}")
                return False
            shutil.copy(str(contigs), str(out_contigs))
            return True
        finally:
            shutil.rmtree(str(work), ignore_errors=True)

    def align_to_insert(self, contigs: Path, out_bam: Path) -> bool:
        """contigs比回insert|Align contigs back to insert"""
        return self.align_to_reference(self.cfg.insert_fasta, contigs, out_bam)

    def align_to_reference(self, ref_fasta, contigs: Path, out_bam: Path) -> bool:
        """contigs比对任意参考|Align contigs to any reference"""
        out_bam.parent.mkdir(parents=True, exist_ok=True)
        ok = self.runner.run_pipeline(
            [[self.cfg.bwa_path, "mem", "-t", str(self.cfg.threads),
              str(ref_fasta), str(contigs)],
             [self.cfg.samtools_path, "sort", "-@", str(self.cfg.threads),
              "-o", str(out_bam), "-"]],
            description=f"contigs比回{Path(ref_fasta).name}|Align contigs to "
                        f"{Path(ref_fasta).name}")
        if not ok:
            return False
        return self.runner.run(
            [self.cfg.samtools_path, "index", str(out_bam)],
            description="索引contigs比对|Index contigs bam")

    def _sam_lines(self, bam: Path) -> List[str]:
        """读bam的SAM行|Read SAM lines from bam"""
        out = self.runner.run_capture(
            [self.cfg.samtools_path, "view", str(bam)],
            description="读取contigs比对|Read contigs alignment")
        return out.splitlines() if out else []

    def _recruit_hits(self, bait_fa: Path, r1: Path, r2: Path) -> List[str]:
        """用诱饵从原始fastq招募read名(去重排序)|Recruit read names via bait (unique, sorted)"""
        out = self.runner.run_pipeline_capture(
            [[self.cfg.bwa_path, "mem", "-t", str(self.cfg.threads),
              str(bait_fa), str(r1), str(r2)],
             [self.cfg.samtools_path, "view", "-q", str(self.cfg.mapq_min),
              "-F", "4", "-"],
             ["cut", "-f1"],
             ["sort", "-u"]],
            description="诱饵招募reads|Bait recruitment")
        if out is None:
            return []
        return [ln for ln in out.splitlines() if ln]

    def _fish_mates(self, names: List[str], mate: int, fastq: Path,
                    pattern_file: Path, out_fastq: Path, append: bool) -> bool:
        """按read名从原始fastq捞mate(委托模块函数)|Fish mates (delegates)"""
        return fish_mates(self.cfg.seqkit_path, self.runner, names, mate,
                          fastq, pattern_file, out_fastq, append)

    def _read_fastq_names(self, paths: List[Path]) -> List[str]:
        """fastq裸read名(去/1/2)|Bare read names from fastq (strip /1,/2)"""
        names = set()
        for p in paths:
            text = p.read_text()
            lines = text.splitlines()
            for i in range(1, len(lines), 4):
                name = lines[i][1:].split()[0]
                if name.endswith("/1") or name.endswith("/2"):
                    name = name[:-2]
                names.add(name)
        return sorted(names)

    # ---------- 主循环|Main loop ----------

    def run(self, sample: str, r1: Path, r2: Path, sc_fq: Path,
            pe1: Path, pe2: Path, walk_dir: Path) -> WalkResult:
        """跑完整步移|Run full walking"""
        walk_dir.mkdir(parents=True, exist_ok=True)
        summary = walk_dir / f"{sample}.walk_summary.tsv"
        master = walk_dir / "master.txt"
        recruited_r1 = walk_dir / "recruited_R1.fastq"
        recruited_r2 = walk_dir / "recruited_R2.fastq"

        # 断点续跑:summary已有记录则从其后继续|Resume from last recorded round
        done_rounds = {}
        if summary.exists():
            done_rounds = {row["round"]: row for row in parse_summary_tsv(summary)}

        if 0 not in done_rounds:
            round0 = walk_dir / "round_0"
            if not self.run_spades(sc_fq, pe1, pe2, round0 / "contigs.fasta"):
                raise RuntimeError("round 0组装失败|Round 0 assembly failed")
            if not self.align_to_insert(round0 / "contigs.fasta",
                                        round0 / "contigs_vs_insert.bam"):
                raise RuntimeError("round 0比对失败|Round 0 alignment failed")
            names = self._read_fastq_names([sc_fq, pe1, pe2])
            master.write_text("\n".join(names) + "\n")
            shutil.copy(str(pe1), str(recruited_r1))
            shutil.copy(str(pe2), str(recruited_r2))
            total, longest, n = contig_stats(round0 / "contigs.fasta")
            self._append_summary(summary, 0, 0, len(names), n, longest, total)
            self.logger.info(
                f"[round 0] contigs={n} longest={longest}bp total={total}bp")
            done_rounds[0] = {"total_bp": total}

        prev_total = contig_stats(_last_round_contigs(walk_dir))[0]

        r = 1
        stop_reason = "max_rounds"
        while r <= self.cfg.max_rounds:
            if r in done_rounds:
                # 已完成轮,跳过|Completed round, skip
                prev_total = contig_stats(walk_dir / f"round_{r}" / "contigs.fasta")[0]
                r += 1
                continue
            curr_contigs = walk_dir / f"round_{r - 1}" / "contigs.fasta"
            curr_bam = walk_dir / f"round_{r - 1}" / "contigs_vs_insert.bam"

            # 1.抽诱饵|Extract bait
            bait_fa = walk_dir / f"bait_{r}.fasta"
            curr_sam = self._sam_lines(curr_bam)
            baits = extract_flank_bait(
                curr_sam, self.cfg.min_softclip, self.cfg.min_unmapped)
            if not baits:
                stop_reason = "no_bait"
                self.logger.info(f"[round {r}] {STOP_REASONS[stop_reason]}")
                break
            with open(bait_fa, "w") as fh:
                for name, seq in baits:
                    fh.write(f">{name}\n{seq}\n")

            # 2.招募|Recruit
            if not self.runner.run([self.cfg.bwa_path, "index", str(bait_fa)],
                                   description=f"诱饵建索引|Index bait round {r}"):
                raise RuntimeError(f"round {r}诱饵索引失败|Bait index failed")
            hits = self._recruit_hits(bait_fa, r1, r2)
            history = set(master.read_text().split())
            new_names = sorted(set(hits) - history)

            if not new_names:
                stop_reason = "no_new_reads"
                self.logger.info(f"[round {r}] {STOP_REASONS[stop_reason]}")
                break
            if len(new_names) > self.cfg.repeat_cap:
                stop_reason = "repeat_cap"
                self.logger.warning(
                    f"[round {r}] 新招{format_number(len(new_names))}条"
                    f"(>{self.cfg.repeat_cap}), {STOP_REASONS[stop_reason]}")
                break

            # 3.并入历史+钓mate|Merge history & fish mates
            history.update(new_names)
            master.write_text("\n".join(sorted(history)) + "\n")
            pat1 = walk_dir / f"pat_r1_{r}.txt"
            pat2 = walk_dir / f"pat_r2_{r}.txt"
            if not self._fish_mates(new_names, 1, r1, pat1, recruited_r1, True):
                raise RuntimeError(f"round {r}钓mate/1失败|Fish mate/1 failed")
            if not self._fish_mates(new_names, 2, r2, pat2, recruited_r2, True):
                raise RuntimeError(f"round {r}钓mate/2失败|Fish mate/2 failed")

            # 4.重组装|Reassemble
            round_dir = walk_dir / f"round_{r}"
            if not self.run_spades(sc_fq, recruited_r1, recruited_r2,
                                   round_dir / "contigs.fasta"):
                raise RuntimeError(f"round {r}组装失败|Assembly failed")
            if not self.align_to_insert(round_dir / "contigs.fasta",
                                        round_dir / "contigs_vs_insert.bam"):
                raise RuntimeError(f"round {r}比对失败|Alignment failed")
            total, longest, n = contig_stats(round_dir / "contigs.fasta")
            self._append_summary(
                summary, r, len(new_names), len(history), n, longest, total)
            self.logger.info(
                f"[round {r}] new={len(new_names)} cum={len(history)} "
                f"contigs={n} longest={longest}bp total={total}bp")

            growth = total - prev_total
            prev_total = total
            # 侧翼延伸判定(向target_flank目标走)|Flank extent (walk toward target)
            curr_flank = max_flank_extent(
                self._sam_lines(round_dir / "contigs_vs_insert.bam"),
                read_fasta_dict(round_dir / "contigs.fasta"))
            reason = decide_stop(len(new_names), growth, True, self.cfg,
                                 current_flank=curr_flank)
            if reason:
                stop_reason = reason
                self.logger.info(
                    f"[round {r}] 侧翼{curr_flank}bp/目标{self.cfg.target_flank}bp, "
                    f"增量{growth}bp, {STOP_REASONS[reason]}")
                break
            self.logger.info(
                f"[round {r}] 侧翼{curr_flank}bp/目标{self.cfg.target_flank}bp,继续步移"
                f"|Flank {curr_flank}/{self.cfg.target_flank}bp, keep walking")
            r += 1
        else:
            self.logger.info(STOP_REASONS["max_rounds"])

        final_contigs = _last_round_contigs(walk_dir)
        total, longest, _ = contig_stats(final_contigs)
        # 完成旗标:summary在round 0后就存在,不能当阶段完成标志|
        # Done flag: summary exists right after round 0, not a stage marker
        (walk_dir / "walk_done.flag").write_text(
            f"stop_reason={stop_reason}\ntotal_bp={total}\nlongest_bp={longest}\n")
        return WalkResult(
            final_contigs=final_contigs, rounds_done=max(done_rounds.keys(),
                                                         default=0) if stop_reason else r,
            stop_reason=stop_reason, total_bp=total, longest_bp=longest)

    @staticmethod
    def _append_summary(summary: Path, rnd: int, new: int, cum: int,
                        n_contigs: int, longest: int, total: int):
        """追加一轮汇总行|Append one summary row"""
        if not summary.exists():
            summary.write_text(
                "round\tnew_reads\tcum_reads\tcontigs\tlongest_bp\ttotal_bp\n")
        with open(summary, "a") as fh:
            fh.write(f"{rnd}\t{new}\t{cum}\t{n_contigs}\t{longest}\t{total}\n")


def _last_round_contigs(walk_dir: Path) -> Path:
    """找最高轮次的contigs|Locate highest-round contigs"""
    rounds = sorted(
        int(p.name.split("_")[1]) for p in walk_dir.glob("round_*")
        if p.is_dir())
    return walk_dir / f"round_{rounds[-1] if rounds else 0}" / "contigs.fasta"
