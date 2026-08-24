"""genome2sv 流程编排器|genome2sv pipeline orchestrator"""
import glob
import os
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple

from .utils import (FaidxReader, ModuleLogger, build_conda_command,
                    build_survivor_input, extract_sv_sequence, format_pav_binary,
                    format_pav_matrix, format_sv_summary_tsv, get_conda_env,
                    gt_present, parse_info_str, parse_svtype_stats, stable_sv_id)


class Genome2SVPipeline:
    """assembly-to-assembly SV calling 编排器|Orchestrator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    # ---------- 命令执行|command execution ----------

    def _run(self, cmd: list, description: str = "") -> bool:
        """执行命令,失败返回 False(不抛)|Run command; return False on failure."""
        if description:
            self.logger.info(f"执行|Executing: {description}")
        self.logger.info(f"命令|Command: {' '.join(str(c) for c in cmd)}")
        try:
            result = subprocess.run(cmd, shell=False, capture_output=True, text=True)
        except FileNotFoundError as e:
            self.logger.error(f"命令未找到|Command not found: {e}")
            return False
        if result.returncode != 0:
            self.logger.error(
                f"命令失败|Command failed (rc={result.returncode}): "
                f"{(result.stderr or '').strip()[:500]}")
            return False
        return True

    # ---------- 步骤 0:参考预处理|reference prep ----------

    def prepare_reference(self) -> bool:
        """软链参考到 reference/ 并 faidx(若缺)|Symlink ref & faidx if missing."""
        ref = self.config.reference_fasta
        src = self.config.reference_path
        if not src:
            self.logger.error("参考路径未解析|Reference path unresolved")
            return False
        try:
            if not os.path.lexists(ref):
                os.symlink(os.path.abspath(src), ref)
                self.logger.info(f"创建参考软链|Symlinked reference: {ref}")
        except OSError as e:
            self.logger.error(f"参考软链失败|Symlink failed: {e}")
            return False
        if Path(self.config.reference_fai).exists():
            self.logger.info("参考 fai 已存在,跳过 faidx|Reference fai exists, skip faidx")
            return True
        cmd = build_conda_command(self.config.samtools_path, ["faidx", ref])
        return self._run(cmd, "samtools faidx 索引参考|samtools faidx reference")

    # ---------- 步骤 1:查询样本(fof 已解析)|queries ----------

    def discover_samples(self) -> List[Tuple[str, str]]:
        """返回查询样本(fof 中除参考外)|Return query samples (fof minus reference)."""
        queries = [(n, p) for n, p in self.config.samples.items()
                   if n != self.config.ref_sample]
        self.logger.info(f"参考样本|Reference: {self.config.ref_sample}")
        self.logger.info(f"查询样本 {len(queries)} 个|{len(queries)} query samples")
        for name, _ in queries:
            self.logger.info(f"  样本|sample: {name}")
        return queries

    # ---------- 步骤 2:比对|alignment ----------

    def _align_command(self, sample: str, query_fasta: str) -> list:
        """构建单 conda run 包裹的 minimap2|samtools 管道
        |Build single-conda-run wrapped minimap2|samtools pipe."""
        bam = str(self.config.alignment_dir / f"{sample}.sorted.bam")
        ref = self.config.reference_fasta
        t = self.config.threads
        script = (
            f"set -eo pipefail; "
            f"minimap2 -ax {self.config.preset} -t {t} {ref} {query_fasta} | "
            f"samtools sort -@ {t} -o {bam} - && "
            f"samtools index {bam}"
        )
        sam_env = get_conda_env(self.config.samtools_path)
        mm_env = get_conda_env(self.config.minimap2_path)
        if sam_env and mm_env and sam_env != mm_env:
            self.logger.warning(
                f"minimap2({mm_env}) 与 samtools({sam_env}) 环境不一致,以 samtools 为准|"
                f"env mismatch minimap2={mm_env} samtools={sam_env}; using samtools env")
        env = sam_env or mm_env
        if env:
            return ["conda", "run", "-n", env, "--no-capture-output",
                    "bash", "-c", script]
        return ["bash", "-c", script]

    def align_sample(self, sample: str, query_fasta: str) -> bool:
        """比对单样本(断点续传)|Align one sample (checkpointed)."""
        bam = self.config.alignment_dir / f"{sample}.sorted.bam"
        bai = self.config.alignment_dir / f"{sample}.sorted.bam.bai"
        if bam.exists() and bai.exists():
            self.logger.info(f"跳过已完成比对|Skipping completed alignment: {sample}")
            return True
        cmd = self._align_command(sample, query_fasta)
        return self._run(cmd, f"minimap2 比对|minimap2 alignment: {sample}")

    # ---------- 步骤 3:SV 调用|SV calling ----------

    @staticmethod
    def _find_svim_vcf(out_dir: str) -> Optional[str]:
        """定位 svim-asm 输出 VCF(glob 兜底)|Locate svim-asm output VCF."""
        hits = sorted(glob.glob(os.path.join(out_dir, "*.svim.vcf")))
        if hits:
            return hits[0]
        hits = sorted(glob.glob(os.path.join(out_dir, "*.vcf")))
        return hits[0] if hits else None

    def call_sv(self, sample: str) -> Optional[str]:
        """svim-asm 调用单样本(断点续传)|Call SVs for one sample (checkpointed).

        Returns:
            VCF 路径或 None(失败)|VCF path or None on failure
        """
        bam = str(self.config.alignment_dir / f"{sample}.sorted.bam")
        out_dir = str(self.config.svim_dir / sample)
        os.makedirs(out_dir, exist_ok=True)
        existing = self._find_svim_vcf(out_dir)
        if existing and os.path.exists(existing):
            self.logger.info(f"跳过已完成 SV 调用|Skipping completed SV calling: {sample}")
            return existing
        args = [self.config.svim_mode, out_dir, bam, self.config.reference_fasta,
                "--min_sv_size", str(self.config.svim_min_sv_size),
                "--sample", sample]
        cmd = build_conda_command(self.config.svim_asm_path, args)
        ok = self._run(cmd, f"svim-asm SV 调用|svim-asm SV calling: {sample}")
        if not ok:
            return None
        vcf = self._find_svim_vcf(out_dir)
        if not vcf:
            self.logger.error(f"svim-asm 未产出 VCF|svim-asm produced no VCF: {sample}")
            return None
        return vcf

    # ---------- 步骤 4:合并 + 统计|merge + stats ----------

    def merge_and_stats(self, sample_vcf_map: dict) -> bool:
        """SURVIVOR 合并 + 统计表|SURVIVOR merge + summary table.

        Args:
            sample_vcf_map: {sample: vcf_path} 仅成功样本|successful samples only
        """
        input_txt = str(self.config.merged_dir / "survivor_input.txt")
        paths = build_survivor_input(sample_vcf_map, input_txt)
        merged_vcf = str(self.config.merged_dir / "pan_sv.survivor.vcf")
        merge_ok = False
        if os.path.exists(merged_vcf):
            # 断点续传:合并是确定性输出,已有结果直接复用
            # |Checkpoint: merge is deterministic; reuse existing result
            self.logger.info(
                "跳过已完成 SURVIVOR 合并|Skipping completed SURVIVOR merge")
            merge_ok = True
        elif not paths:
            self.logger.warning(
                "无可用样本 VCF,跳过 SURVIVOR 合并|No sample VCFs; skip SURVIVOR merge")
        else:
            # SURVIVOR v1.0.7: merge filelist max_dist min_support type strand est_dist min_sv_len out
            args = ["merge", input_txt, str(self.config.max_dist),
                    str(self.config.min_support), str(self.config.survivor_type),
                    str(self.config.survivor_strand), str(self.config.est_dist),
                    str(self.config.min_sv_length), merged_vcf]
            cmd = build_conda_command(self.config.survivor_path, args)
            merge_ok = self._run(cmd, "SURVIVOR 合并|SURVIVOR merge")
            if merge_ok and os.path.exists(merged_vcf):
                self._bcftools_stats(merged_vcf)
            elif not merge_ok:
                self.logger.error("SURVIVOR 合并失败|SURVIVOR merge failed")

        # 统计表(始终生成)|summary always generated
        rows = [(s, parse_svtype_stats(sample_vcf_map[s]))
                for s in sorted(sample_vcf_map)]
        if merge_ok and os.path.exists(merged_vcf):
            rows.append(("merged", parse_svtype_stats(merged_vcf)))
        tsv = format_sv_summary_tsv(rows)
        summary_path = self.config.stats_dir / "sv_summary.tsv"
        summary_path.write_text(tsv)
        self.logger.info(f"统计表已生成|Summary written: {summary_path}")
        return True

    def _bcftools_stats(self, merged_vcf: str) -> None:
        """bcftools stats 写入文件|Write bcftools stats to file."""
        cmd = build_conda_command(self.config.bcftools_path, ["stats", merged_vcf])
        self.logger.info(f"命令|Command: {' '.join(str(c) for c in cmd)}")
        try:
            r = subprocess.run(cmd, shell=False, capture_output=True, text=True)
            (self.config.stats_dir / "pan_sv.survivor.stats").write_text(r.stdout)
        except FileNotFoundError as e:
            self.logger.warning(f"bcftools stats 跳过|bcftools stats skipped: {e}")

    # ---------- 步骤 5:SV 序列 + PAV 矩阵|sequences + PAV matrix ----------

    def generate_downstream_outputs(self, merged_vcf: str) -> bool:
        """从 merged VCF 生成 SV 序列 FASTA + PAV 矩阵|Sequences FASTA + PAV from merged VCF.

        SV 序列:INS 取 ALT(剥 anchor)、DEL 取 REF、INV 取 ALT(已 revcomp)、
        DUP 按参考 [POS,END] 提取重复单元;BND/无法提取者跳过并计数。
        PAV:GT 含 1 记 1,./. 记 0;与序列 FASTA 共用自增 sv_id 便于交叉引用。
        |Sequences: INS from ALT (anchor stripped), DEL from REF, INV from ALT
        (already revcomp), DUP as reference [POS,END] unit; BND/unextractable
        skipped with counts. PAV: GT with allele 1 → 1 else 0; shares the
        auto-increment sv_id with the FASTA for cross-reference.
        """
        seq_fa = self.config.sv_seq_dir / "pan_sv.sequences.fa"
        pav_tsv = self.config.stats_dir / "pav_matrix.tsv"
        pav_bin = self.config.stats_dir / "pav_binary.tsv"
        if all(p.exists() for p in (seq_fa, pav_tsv, pav_bin)):
            self.logger.info(
                "跳过已完成下游输出|Skipping completed downstream outputs "
                "(SV sequences + PAV)")
            return True
        if not os.path.exists(merged_vcf):
            self.logger.warning(
                f"merged VCF 不存在,跳过 SV 序列与 PAV 输出|merged VCF missing; "
                f"skip SV sequences & PAV: {merged_vcf}")
            return True

        reader = None
        try:
            reader = FaidxReader(self.config.reference_fasta)
        except (FileNotFoundError, OSError) as e:
            self.logger.warning(
                f"参考 fai 不可用,坐标提取(DUP/符号化回退)将跳过|Reference fai "
                f"unavailable; region extraction (DUP/symbolic fallback) skipped: {e}")

        sample_names: List[str] = []
        type_counters: dict = {}
        pav_rows: List[tuple] = []
        fa_records: List[Tuple[str, str]] = []
        skipped_seq = 0
        with open(merged_vcf) as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    sample_names = line.rstrip("\n").split("\t")[9:]
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue
                chrom, pos = fields[0], int(fields[1])
                ref, alt = fields[3], fields[4]
                info = parse_info_str(fields[7])
                svtype = info.get("SVTYPE", "UNKNOWN")
                end = info.get("END", pos)
                svlen = info.get("SVLEN", ".")
                n = type_counters.get(svtype, 0) + 1
                type_counters[svtype] = n
                sv_id = stable_sv_id(svtype, n)
                pav = [gt_present(f) for f in fields[9:]] if len(fields) > 9 else []
                pav_rows.append((sv_id, chrom, pos, end, svtype, svlen, pav))
                result = extract_sv_sequence(svtype, chrom, pos, ref, alt,
                                             info, reader)
                if result is None:
                    skipped_seq += 1
                    continue
                seq, source = result
                support = ",".join(s for s, p in zip(sample_names, pav) if p)
                header = (f">{sv_id} type={svtype.split(':', 1)[0]} chrom={chrom} "
                          f"pos={pos} end={end} len={len(seq)} source={source} "
                          f"samples={support}")
                fa_records.append((header, seq))
        if reader is not None:
            reader.close()

        # FASTA(60 列换行)|FASTA wrapped at 60 columns
        with open(seq_fa, "w") as out:
            for header, seq in fa_records:
                out.write(header + "\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + "\n")
        pav_tsv.write_text(format_pav_matrix(pav_rows, sample_names))
        pav_bin.write_text(format_pav_binary(pav_rows, sample_names))

        self.logger.info(
            f"SV 序列已生成|SV sequences written: {seq_fa} "
            f"(提取|extracted {len(fa_records)}, 跳过|skipped {skipped_seq})")
        for svtype in sorted(type_counters):
            self.logger.info(
                f"  {svtype}: {type_counters[svtype]} 条(进 PAV)|records (in PAV)")
        self.logger.info(f"PAV 矩阵已生成|PAV matrices written: {pav_tsv}, {pav_bin}")
        return True

    # ---------- 主流程|main run ----------

    def run(self) -> int:
        """端到端运行|Run end-to-end. Returns exit code 0/1."""
        from datetime import datetime
        start_time = datetime.now()
        if not self.prepare_reference():
            self.logger.error("参考预处理失败,终止|Reference prep failed; abort")
            return 1
        queries = self.discover_samples()

        sample_vcf_map = {}
        for sample, query_fasta in queries:
            if not self.align_sample(sample, query_fasta):
                self.logger.error(f"样本比对失败,跳过|Alignment failed; skip: {sample}")
                continue
            vcf = self.call_sv(sample)
            if not vcf:
                self.logger.error(f"样本 SV 调用失败,跳过|SV calling failed; skip: {sample}")
                continue
            sample_vcf_map[sample] = vcf

        if not sample_vcf_map:
            self.logger.error("全部样本失败,无合并输入|All samples failed; nothing to merge")
            return 1

        self.merge_and_stats(sample_vcf_map)
        merged_vcf = str(self.config.merged_dir / "pan_sv.survivor.vcf")
        if os.path.exists(merged_vcf):
            self.generate_downstream_outputs(merged_vcf)
        else:
            self.logger.warning(
                "无 merged VCF,跳过 SV 序列与 PAV 输出|No merged VCF; "
                "skip SV sequences & PAV output")
        from .utils import write_software_versions
        write_software_versions(
            self.config, self.logger,
            str(self.config.info_dir / "software_versions.yml"),
            start_time=start_time)
        self.logger.info("genome2sv 流程完成|genome2sv pipeline completed")
        return 0
