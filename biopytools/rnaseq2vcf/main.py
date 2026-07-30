"""rnaseq2vcf 主入口|orchestrator + argparse entry (per-sample gVCF → joint calling → one multi-sample VCF)

流程:HISAT2 比对 → 每样本 gVCF(HaplotypeCaller -ERC GVCF) → CombineGVCFs → GenotypeGVCFs
→ VariantFiltration → bcftools PASS → 一个多样本 VCF。注释由用户手动跑 biopytools annovar。|
Per-sample gVCF → joint genotyping → one multi-sample VCF. Annotate via biopytools annovar.
"""

import argparse
import os
import subprocess
import sys

from .config import Rnaseq2vcfConfig
from .utils import (Rnaseq2vcfLogger, CommandRunner, CheckpointManager,
                    SystemChecker, build_conda_command, discover_samples)
from .data_processing import GenomeIndexer, QualityController, Aligner, Caller, JointCaller


class Rnaseq2vcfProcessor:
    """端到端编排器(多样本联合 calling)|End-to-end orchestrator (multi-sample joint calling)"""

    def __init__(self, config: Rnaseq2vcfConfig):
        self.config = config
        self.logger_manager = Rnaseq2vcfLogger(
            os.path.join(config.output_dir, "99_logs"),
            log_file=config.log_file, log_level=config.log_level)
        self.logger = self.logger_manager.get_logger()
        self.runner = CommandRunner(self.logger, config.output_dir, dry_run=config.dry_run)
        self.checkpoints = CheckpointManager(
            os.path.join(config.output_dir, "00_pipeline_info", "checkpoints"), self.logger)

        self.indexer = GenomeIndexer(config, self.logger, self.runner)
        self.qc_controller = QualityController(config, self.logger, self.runner)
        self.aligner = Aligner(config, self.logger, self.runner)
        self.caller = Caller(config, self.logger, self.runner)
        self.joint_caller = JointCaller(config, self.logger, self.runner)

    def _run_shared_index(self) -> bool:
        # force=True 时绕过 checkpoint,强制重建索引(否则 -f 对共享索引无效)|
        # force bypasses the checkpoint so -f actually rebuilds the shared index
        if (self.config.enable_checkpoint and not self.config.force
                and self.checkpoints.exists("genome_index")):
            self.logger.info("跳过已完成步骤|Skipping completed step: genome_index")
            return True
        ok = self.indexer.run()
        if ok and self.config.enable_checkpoint:
            self.checkpoints.create("genome_index")
        return ok

    # ---- 单样本包装|per-sample wrappers(供测试 mock)----
    def qc(self, sample, r1, r2):
        if self.config.clean_fastq_dir or self.config.skip_qc:
            return r1, r2
        return self.qc_controller.run_sample(sample, r1, r2)

    def align(self, sample, r1, r2):
        return self.aligner.run_sample(sample, r1, r2)

    def call(self, sample, bam):
        return self.caller.run_sample(sample, bam)

    def joint(self, gvcfs):
        return self.joint_caller.run(gvcfs)

    def _run_one_sample(self, sample, r1, r2):
        """单样本走 QC→align→call(gVCF),返回 gVCF 路径|per-sample QC→align→gVCF, returns gVCF path"""
        self.logger.info(f"==== 处理样本|Processing sample: {sample} ====")
        r1c, r2c = self.qc(sample, r1, r2)
        bam = self.align(sample, r1c, r2c)
        gvcf = self.call(sample, bam)
        self.logger.info(f"样本 gVCF 完成|Sample gVCF done: {sample} → {gvcf}")
        return gvcf

    def _pre_flight(self) -> bool:
        cfg = self.config
        for desc, t in (("hisat2", cfg.hisat2_path), ("gatk", cfg.gatk_path),
                        ("samtools", cfg.samtools_path), ("bcftools", cfg.bcftools_path)):
            if not SystemChecker.check_command_exists(t, self.logger):
                return False
        SystemChecker.check_disk_space(cfg.output_dir, 200, self.logger)
        return True

    def run(self):
        cfg = self.config
        self.logger.info("rnaseq2vcf 启动|Starting rnaseq2vcf pipeline (joint calling)")
        if not self._pre_flight():
            self.logger.error("预检失败,退出|Pre-flight failed, exiting")
            return False
        try:
            cfg.validate()
        except ValueError as e:
            self.logger.error(str(e))
            return False

        if not self._run_shared_index():
            return False
        if cfg.step == 0:
            self.logger.info("仅构建共享索引,退出|Index-only requested, exiting")
            return True

        input_dir = cfg.clean_fastq_dir or cfg.raw_fastq_dir
        samples = discover_samples(input_dir, cfg.read1_pattern, cfg.read2_pattern)
        self.logger.info(f"发现样本|Found {len(samples)} samples: {[s[0] for s in samples]}")

        failed = []
        sample_gvcfs = []
        for sample, r1, r2 in samples:
            try:
                gvcf = self._run_one_sample(sample, r1, r2)
                sample_gvcfs.append(gvcf)
            except Exception as e:
                self.logger.error(f"样本失败,继续下一个|Sample failed, continuing: {sample}: {e}")
                failed.append(sample)

        # 联合 calling:所有样本 gVCF → 一个多样本 VCF|joint calling: all gVCFs → one multi-sample VCF
        joint_vcf = None
        if sample_gvcfs:
            try:
                joint_vcf = self.joint(sample_gvcfs)
                self.logger.info(f"联合 VCF 完成|Joint VCF done: {joint_vcf}")
            except Exception as e:
                self.logger.error(f"联合 calling 失败|Joint calling failed: {e}")
                failed.append("joint_calling")

        self._write_report(samples, failed, joint_vcf)
        if failed:
            self.logger.warning(f"失败项|Failed: {failed}")
        return len(failed) == 0

    def _count_records(self, vcf: str):
        """统计 VCF 变异记录数(流式,内存安全)|count variant records (streaming, memory-safe)"""
        cmd = build_conda_command(self.config.bcftools_path, ['view', '-H', vcf])
        self.logger.info(f"统计变异数|Counting variants: {vcf}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            n = sum(1 for _ in proc.stdout)
            proc.wait()
            if proc.returncode != 0:
                # bcftools 失败(VCF 损坏等)时 n 可能是部分计数,不可信|partial count on failure is unreliable
                self.logger.warning(f"统计命令失败|Count command failed (rc={proc.returncode}): {vcf}")
                return None
            return n
        except Exception as e:
            self.logger.warning(f"统计失败|Count failed ({vcf}): {e}")
            return None

    def _per_sample_counts(self, vcf: str, sample_names):
        """每样本非参考基因型变异数(bcftools stats PSC)|per-sample non-ref variant count"""
        cmd = build_conda_command(self.config.bcftools_path, ['stats', '-s', '-', vcf])
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        counts = {s: 0 for s in sample_names}
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
            for line in r.stdout.splitlines():
                if line.startswith('PSC'):
                    f = line.split('\t')
                    if len(f) >= 6:
                        # PSC: id sample nRefHom nNonRefHom nHet ...
                        counts[f[2]] = int(f[4]) + int(f[5])
        except Exception as e:
            self.logger.warning(f"每样本统计失败|Per-sample stats failed ({vcf}): {e}")
        return counts

    def _write_report(self, samples, failed, joint_vcf):
        cfg = self.config
        report = os.path.join(cfg.output_dir, "ANALYSIS_REPORT.txt")
        joint_dir = os.path.join(cfg.output_dir, "04_joint")
        pre_vcf = os.path.join(joint_dir, "joint.vcf.gz")
        filt_vcf = os.path.join(joint_dir, "joint.filtered.vcf.gz")
        pass_vcf = joint_vcf or os.path.join(joint_dir, "all_samples.pass.vcf.gz")

        n_pre = self._count_records(pre_vcf) if os.path.exists(pre_vcf) else None
        n_pass = self._count_records(pass_vcf) if os.path.exists(pass_vcf) else None
        per_sample = (self._per_sample_counts(pass_vcf, [s[0] for s in samples])
                      if os.path.exists(pass_vcf) else {})

        lines = []
        lines.append("rnaseq2vcf 分析报告|Analysis Report (joint calling; ANNOVAR 由用户手动运行|run annovar manually)")
        lines.append("=" * 64)
        lines.append(f"参考基因组|Reference genome: {cfg.ref_genome_fa}")
        lines.append(f"GFF3(剪接位点)|GFF3 (splice sites): {cfg.gff3_file or '(未提供|none)'}")
        lines.append(f"样本数|Sample count: {len(samples)}")
        lines.append(f"样本|Samples: {' '.join(s[0] for s in samples)}")
        lines.append(f"失败项|Failed: {failed}")
        lines.append("")
        lines.append("变异过滤统计|Variant filtering (joint VCF; FS>30 / QD<2 / cluster 3 @ 35bp)")
        lines.append("-" * 64)
        if n_pre is not None and n_pass is not None:
            out = n_pre - n_pass
            pct = (out / n_pre * 100) if n_pre else 0.0
            lines.append(f"过滤前变异总数|Total before filter  (joint.vcf.gz):           {n_pre}")
            lines.append(f"过滤后保留 PASS|Retained PASS (all_samples.pass.vcf.gz):     {n_pass}")
            lines.append(f"被过滤掉|Filtered out:                                   {out} ({pct:.1f}%)")
        else:
            lines.append("变异统计不可用(联合 VCF 未生成)|Counts unavailable (joint VCF missing)")
        lines.append("")
        if per_sample:
            lines.append("每样本 PASS 变异数(非参考基因型)|Per-sample PASS variants (non-ref GT):")
            for s, _, _ in samples:
                lines.append(f"  {s}: {per_sample.get(s, 0)}")
            lines.append("")
        lines.append("下游分析文件|Downstream analysis file")
        lines.append("=" * 64)
        if pass_vcf and os.path.exists(pass_vcf):
            lines.append(f"  ★ {pass_vcf}")
            lines.append("    (含全部样本,已过滤;送 annovar / 其他下游|all samples, filtered; for annovar/downstream)")
            lines.append(f"    索引|index: {pass_vcf}.tbi")
            lines.append(f"    复查被滤变异|review filtered: {filt_vcf} (FILTER 列含 FS/QD/SNPCluster)")
        else:
            lines.append("  [最终 VCF 缺失|final VCF missing]")
        with open(report, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines) + '\n')
        self.logger.info(f"报告已写|Report written: {report}")


def parse_arguments(argv=None):
    p = argparse.ArgumentParser(
        description="转录组变异检测(多样本联合 calling)|RNA-seq variant calling (multi-sample joint calling)")
    p.add_argument('-g', '--genome', required=True, help='参考基因组 FASTA|Reference genome FASTA')
    p.add_argument('--gff3', help='参考 GFF3(可选,HISAT2 剪接位点)|Reference GFF3 (optional, HISAT2 splice sites)')
    p.add_argument('-i', '--input', help='原始 FASTQ 目录|Raw FASTQ dir (runs fastp)')
    p.add_argument('--clean-fastq-dir', help='已清洗 FASTQ 目录(跳过 QC)|Clean FASTQ dir (skip QC)')
    p.add_argument('-o', '--output-dir', default='.', help='输出目录|Output dir')
    p.add_argument('-t', '--threads', type=int, default=12, help='线程数|Threads (default 12)')
    p.add_argument('--min-conf', type=int, default=20, help='HaplotypeCaller 最小置信度|min confidence')
    p.add_argument('--fs-threshold', type=float, default=30.0)
    p.add_argument('--qd-threshold', type=float, default=2.0)
    p.add_argument('--cluster-window', type=int, default=35)
    p.add_argument('--cluster-size', type=int, default=3)
    p.add_argument('--read1-pattern', help='R1 后缀(默认自动识别 _1.clean.fq.gz/_1.fq.gz 等)|R1 suffix (auto)')
    p.add_argument('--read2-pattern', help='R2 后缀(默认自动识别)|R2 suffix (auto)')
    p.add_argument('-s', '--step', type=int, choices=[0, 1, 2, 3, 4],
                   help='0=仅建索引|index only;省略=全流程|omit for full pipeline')
    p.add_argument('--no-checkpoint', action='store_true')
    p.add_argument('-f', '--force', action='store_true')
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--skip-qc', action='store_true')
    p.add_argument('--log-file')
    p.add_argument('--log-level', default='INFO')
    return p.parse_args(argv)


def main():
    args = parse_arguments()
    try:
        cfg = Rnaseq2vcfConfig(
            ref_genome_fa=args.genome, gff3_file=args.gff3, output_dir=args.output_dir,
            raw_fastq_dir=args.input, clean_fastq_dir=args.clean_fastq_dir,
            threads=args.threads, min_conf=args.min_conf,
            fs_threshold=args.fs_threshold, qd_threshold=args.qd_threshold,
            cluster_window=args.cluster_window, cluster_size=args.cluster_size,
            read1_pattern=args.read1_pattern, read2_pattern=args.read2_pattern,
            step=args.step,
            enable_checkpoint=not args.no_checkpoint, dry_run=args.dry_run,
            force=args.force, skip_qc=args.skip_qc, log_file=args.log_file, log_level=args.log_level)
        proc = Rnaseq2vcfProcessor(cfg)
        ok = proc.run()
        sys.exit(0 if ok else 1)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
