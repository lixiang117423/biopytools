"""genome2sv 流程编排器|genome2sv pipeline orchestrator"""
import glob
import os
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple

from .utils import (ModuleLogger, build_conda_command, get_conda_env,
                    build_survivor_input, format_sv_summary_tsv, parse_svtype_stats)


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
        if not paths:
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
        from .utils import write_software_versions
        write_software_versions(
            self.config, self.logger,
            str(self.config.info_dir / "software_versions.yml"),
            start_time=start_time)
        self.logger.info("genome2sv 流程完成|genome2sv pipeline completed")
        return 0
