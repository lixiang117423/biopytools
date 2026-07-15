"""
ps-gene-anno 命令行入口 + 流程编排(GFF3)|CLI entry + pipeline (GFF3)
BRAKER 后效应子查漏补缺|Post-BRAKER effector gap-filling
"""

import argparse
import os
import sys
from typing import List

from .config import PsGeneAnnoConfig
from .utils import PsGeneAnnoLogger
from .evidence import run_miniprot, parse_miniprot_gff3
from .gap_analysis import (
    parse_braker_gff3, detect_gaps, detect_merged_genes, parse_repeat_out,
    dedupe_hits)
from .model_build import qc_filter, build_gene_models
from .merge import merge_gff3


def parse_arguments():
    """解析命令行参数|Parse CLI arguments"""
    parser = argparse.ArgumentParser(
        description="ps-gene-anno: BRAKER 后效应子查漏补缺(GFF3)"
        "|Post-BRAKER effector gap-filling (GFF3)")
    parser.add_argument('-g', '--genome', required=True,
                        help='未mask原始基因组|Unmasked raw genome')
    parser.add_argument('-b', '--braker-gff3', required=True,
                        help='BRAKER输出GFF3|BRAKER output GFF3')
    parser.add_argument('-p', '--prot-seq', required=True,
                        help='近缘蛋白|Protein evidence')
    parser.add_argument('-o', '--output-dir', required=True,
                        help='输出目录|Output directory')
    parser.add_argument('--rnaseq-bam', help='RNA-seq BAM(逗号分隔)|RNA-seq BAMs')
    parser.add_argument('--isoseq-bam', help='三代BAM|Long-read BAM')
    parser.add_argument('--repeat-out', help='RepeatMasker .out|RepeatMasker out')
    parser.add_argument('--prefix', help='输出前缀(默认genome stem)|Output prefix')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Threads')
    # 质控
    parser.add_argument('--gap-min-identity', type=float, default=70)
    parser.add_argument('--gap-min-coverage', type=float, default=80)
    parser.add_argument('--gap-min-cds-len', type=int, default=100)
    parser.add_argument('--overlap-cutoff', type=float, default=0)
    parser.add_argument('--no-require-complete-orf', dest='require_complete_orf',
                        action='store_false', default=True)
    parser.add_argument('--te-overlap-cutoff', type=float, default=50)
    parser.add_argument('--exclude-te-gap', action='store_true',
                        help='质控排除TE区gap(默认不排,疫霉效应子常在TE区)|exclude TE-overlap gaps')
    # 合并拆分
    parser.add_argument('--no-split', dest='enable_split',
                        action='store_false', default=True)
    parser.add_argument('--split-min-hits', type=int, default=2)
    parser.add_argument('--split-min-copy-coverage', type=float, default=80)
    # 步骤控制
    parser.add_argument('--skip-evidence-scan', action='store_true')
    parser.add_argument('--skip-gap-analysis', action='store_true')
    parser.add_argument('--skip-merge', action='store_true')
    return parser.parse_args()


class PsGeneAnnoRunner:
    """ps-gene-anno 流程编排|Pipeline orchestrator"""

    def __init__(self, config: PsGeneAnnoConfig, logger=None):
        self.config = config
        if logger is None:
            log_file = os.path.join(config.log_dir, 'ps_gene_anno.log')
            self.logger = PsGeneAnnoLogger(log_file).get_logger()
        else:
            self.logger = logger
        from .utils import CommandRunner
        self.cmd_runner = CommandRunner(self.logger, config.output_dir)
        # 输出路径(GFF3)|output paths (GFF3)
        self.miniprot_gff = os.path.join(
            config.evidence_dir, f"{config.prefix}.miniprot.gff3")
        self.gap_filled_gff3 = os.path.join(
            config.gap_filled_dir, f"{config.prefix}.gap_filled.gff3")
        self.merged_gff3 = os.path.join(
            config.merged_dir, f"{config.prefix}.merged.gff3")

    def _step_evidence_scan(self):
        """Step1: miniprot 证据扫描(断点续传)|miniprot scan (checkpoint)"""
        if self.config.skip_evidence_scan or os.path.exists(self.miniprot_gff):
            self.logger.info("跳过证据扫描|Skipping evidence scan")
            return
        run_miniprot(self.config.genome, self.config.prot_seq,
                     self.miniprot_gff, self.config,
                     self.cmd_runner, self.logger)

    def _load_hits(self):
        """加载并预过滤 miniprot 命中|Load + pre-filter hits"""
        hits = parse_miniprot_gff3(
            self.miniprot_gff,
            self.config.gap_min_identity, self.config.gap_min_coverage)
        self.logger.info(f"miniprot 命中(过滤后)|hits filtered: {len(hits)}")
        hits = dedupe_hits(hits)   # 全 prot: 同位置多 query 合并|dedup same-locus
        self.logger.info(f"miniprot 命中(去重后)|hits deduped: {len(hits)}")
        return hits

    def run(self):
        """运行完整流程|Run full pipeline"""
        self.logger.info("=" * 70)
        self.logger.info("ps-gene-anno: BRAKER 后效应子查漏补缺(GFF3)|Post-BRAKER gap-filling")
        self.logger.info("=" * 70)

        # Step1 证据扫描|evidence scan
        self._step_evidence_scan()
        if not os.path.exists(self.miniprot_gff):
            self.logger.error("证据扫描未产出|Evidence scan produced nothing")
            return None
        hits = self._load_hits()

        # Step2 漏检/合并分析|gap/merged analysis
        if self.config.skip_gap_analysis:
            self.logger.info("跳过分析|Skipping gap analysis")
            return None
        braker_genes = parse_braker_gff3(self.config.braker_gff3)
        self.logger.info(f"braker 基因|braker genes: {len(braker_genes)}")

        gaps = detect_gaps(hits, braker_genes, self.config.overlap_cutoff)
        self.logger.info(f"漏检候选|gap candidates: {len(gaps)}")

        merged = []
        if self.config.enable_split:
            merged = detect_merged_genes(
                hits, braker_genes,
                self.config.split_min_hits, self.config.split_min_copy_coverage)
            self.logger.info(f"错误合并基因|merged genes: {len(merged)}")

        # 排除已归入 merged 基因的漏检(避免重复补)|exclude hits in merged genes
        merged_bounds = [(g.chrom, g.start, g.end) for g, _ in merged]
        gaps = [h for h in gaps
                if not any(h.chrom == c and h.start >= s and h.end <= e
                           for c, s, e in merged_bounds)]
        merged_hits = [h for _, hs in merged for h in hs]
        all_candidate_hits = gaps + merged_hits

        # Step3 质控 + 建模型(GFF3)|QC + model build
        repeat_regions = parse_repeat_out(self.config.repeat_out) \
            if self.config.repeat_out else {}
        if not self.config.repeat_out:
            self.logger.warning("未提供 repeat_out, 跳过真TE区排除"
                                "|No repeat_out, skipping real-TE exclusion")
        passed = qc_filter(all_candidate_hits, self.config, repeat_regions)
        self.logger.info(f"质控通过|QC passed: {len(passed)}/{len(all_candidate_hits)}")

        gap_lines = build_gene_models(passed, self.config.prefix)
        with open(self.gap_filled_gff3, 'w') as f:
            f.write("\n".join(gap_lines))
            if gap_lines:
                f.write("\n")
        self.logger.info(f"gap_filled 写出|gap_filled written: {self.gap_filled_gff3}")

        # Step4 合并(GFF3)|merge
        if self.config.skip_merge:
            self.logger.info("跳过合并|Skipping merge")
            return self.gap_filled_gff3
        merged_gene_ids = {g.gene_id for g, _ in merged}
        merge_gff3(self.config.braker_gff3, gap_lines,
                   merged_gene_ids, self.merged_gff3)
        self.logger.info(f"merged 写出|merged written: {self.merged_gff3}")

        # Step5 gap 验证报告(蛋白+RNA-seq+TE)|gap evidence report
        report_tsv = os.path.join(self.config.gap_dir,
                                  f"{self.config.prefix}.gap_report.tsv")
        from .report import build_gap_report
        build_gap_report(passed, self.config.prefix, self.config.rnaseq_bam,
                         self.config.repeat_out, report_tsv,
                         self.config, self.cmd_runner, self.logger)

        self.logger.info("=" * 70)
        self.logger.info("ps-gene-anno 完成|ps-gene-anno done")
        self.logger.info("=" * 70)
        return self.merged_gff3


def main():
    """主入口|Main entry"""
    args = parse_arguments()
    rnaseq = args.rnaseq_bam.split(',') if args.rnaseq_bam else None
    try:
        config = PsGeneAnnoConfig(
            genome=args.genome, braker_gff3=args.braker_gff3,
            prot_seq=args.prot_seq, output_dir=args.output_dir,
            rnaseq_bam=rnaseq, isoseq_bam=args.isoseq_bam,
            repeat_out=args.repeat_out, prefix=args.prefix,
            gap_min_identity=args.gap_min_identity,
            gap_min_coverage=args.gap_min_coverage,
            gap_min_cds_len=args.gap_min_cds_len,
            overlap_cutoff=args.overlap_cutoff,
            require_complete_orf=args.require_complete_orf,
            te_overlap_cutoff=args.te_overlap_cutoff,
            exclude_te_gap=args.exclude_te_gap,
            enable_split=args.enable_split,
            split_min_hits=args.split_min_hits,
            split_min_copy_coverage=args.split_min_copy_coverage,
            threads=args.threads,
            skip_evidence_scan=args.skip_evidence_scan,
            skip_gap_analysis=args.skip_gap_analysis,
            skip_merge=args.skip_merge,
        )
        config.validate()
        runner = PsGeneAnnoRunner(config)
        result = runner.run()
        if result:
            sys.exit(0)
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
