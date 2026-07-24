"""rnaseq2vcf 主入口|orchestrator + argparse entry (VCF-only, no annotation)

流程止于生成 PASS VCF;ANNOVAR 注释由用户后续手动调用 `biopytools annovar`。|
Pipeline stops at PASS VCF; ANNOVAR annotation is run manually via `biopytools annovar`.
"""

import argparse
import os
import sys

from .config import Rnaseq2vcfConfig
from .utils import (Rnaseq2vcfLogger, CommandRunner, CheckpointManager,
                    SystemChecker, discover_samples)
from .data_processing import GenomeIndexer, QualityController, Aligner, Caller, VariantFilter


class Rnaseq2vcfProcessor:
    """端到端编排器(到 VCF)|End-to-end orchestrator (up to VCF)"""

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
        self.filterer = VariantFilter(config, self.logger, self.runner)

    def _run_shared_index(self) -> bool:
        if self.config.enable_checkpoint and self.checkpoints.exists("genome_index"):
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

    def filter(self, sample, raw_vcf):
        return self.filterer.run_sample(sample, raw_vcf)

    def _run_one_sample(self, sample, r1, r2):
        self.logger.info(f"==== 处理样本|Processing sample: {sample} ====")
        r1c, r2c = self.qc(sample, r1, r2)
        bam = self.align(sample, r1c, r2c)
        raw_vcf = self.call(sample, bam)
        pass_vcf = self.filter(sample, raw_vcf)
        self.logger.info(f"样本完成|Sample done: {sample} → PASS VCF: {pass_vcf}")

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
        self.logger.info("rnaseq2vcf 启动|Starting rnaseq2vcf pipeline (VCF-only)")
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
        for sample, r1, r2 in samples:
            try:
                self._run_one_sample(sample, r1, r2)
            except Exception as e:
                self.logger.error(f"样本失败,继续下一个|Sample failed, continuing: {sample}: {e}")
                failed.append(sample)

        self._write_report(samples, failed)
        if failed:
            self.logger.warning(f"失败样本|Failed samples: {failed}")
        return len(failed) == 0

    def _write_report(self, samples, failed):
        cfg = self.config
        report = os.path.join(cfg.output_dir, "ANALYSIS_REPORT.txt")
        lines = ["rnaseq2vcf 分析报告|Analysis Report (VCF-only; ANNOVAR 由用户手动运行|ANNOVAR run manually)",
                 f"参考基因组|Reference genome: {cfg.ref_genome_fa}",
                 f"GFF3(剪接位点)|GFF3 (splice sites): {cfg.gff3_file}",
                 f"样本数|Sample count: {len(samples)}",
                 f"失败样本|Failed: {failed}",
                 "", "PASS VCF:"]
        for s, _, _ in samples:
            vcf = os.path.join(cfg.output_dir, s, "04_filter", f"{s}.pass.vcf.gz")
            lines.append(f"  {s}: {vcf} {'[OK]' if os.path.exists(vcf) else '[MISSING]'}")
        with open(report, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines) + '\n')
        self.logger.info(f"报告已写|Report written: {report}")


def parse_arguments(argv=None):
    p = argparse.ArgumentParser(
        description="转录组变异检测(到 VCF)|RNA-seq variant calling (up to VCF)")
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
