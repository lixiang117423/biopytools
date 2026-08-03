"""
annorefine 主程序|annorefine Main Entry
端到端基因组注释: BRAKER 注释 + 同源查缺补漏 → 整合 GFF3
|End-to-end genome annotation: BRAKER + homology gap-filling → integrated GFF3

输入: 未mask基因组 + 转录组(RNA-seq 目录) + 同源蛋白 → 一条命令跑完
|Input: unmasked genome + transcriptome (RNA-seq dirs) + homolog proteins → one command
"""

import argparse
import os
import subprocess
import sys

from .config import AnnorefineConfig
from .utils import AnnorefineLogger
from .evidence import run_miniprot, parse_miniprot_gff3, load_fasta
from .gap_analysis import (
    parse_braker_gff3, detect_gaps, detect_merged_genes, parse_repeat_out,
    dedupe_hits)
from .model_build import qc_filter, qc_filter_small, build_gene_models
from .expression import prepare_unique_bam, compute_hit_depth_breadth
from .merge import merge_gff3


class AnnorefineRunner:
    """查漏补缺引擎(同源补漏 + 合并拆分 + ORF/表达质控)
    |Gap-filling engine: homology gap-fill + merged-gene split + ORF/expression QC"""

    def __init__(self, config: AnnorefineConfig, logger=None):
        self.config = config
        if logger is None:
            log_file = os.path.join(config.log_dir, 'annorefine.log')
            self.logger = AnnorefineLogger(log_file).get_logger()
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

    def _probe_version(self, tool_path: str) -> str:
        """获取工具版本(conda 包装, §13)|probe tool version (conda-wrapped)"""
        if not tool_path:
            return 'not_configured'
        from .utils import get_conda_env
        env = get_conda_env(tool_path)
        if env:
            cmd = ['conda', 'run', '-n', env, '--no-capture-output',
                   tool_path, '--version']
        else:
            cmd = [tool_path, '--version']
        try:
            self.logger.info(f"命令|Command: {' '.join(cmd)}")
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            ver = (r.stdout or r.stderr or '').strip().split('\n')[0]
            return ver or 'unknown'
        except Exception as e:
            self.logger.warning(
                f"获取版本失败|version probe failed ({tool_path}): {e}")
            return 'unknown'

    def _generate_software_versions(self):
        """生成 software_versions.yml 到 00_pipeline_info(§12.5)|write versions yml"""
        import yaml
        try:
            from . import __version__ as mod_ver
        except Exception:
            mod_ver = 'unknown'
        # 查漏补缺阶段用到的工具(miniprot/samtools/stringtie)
        # |tools used in gap-filling phase
        tools = {
            'miniprot': self.config.miniprot_bin,
            'samtools': self.config.samtools_bin,
            'stringtie': self.config.stringtie_bin,
        }
        versions = {name: {'version': self._probe_version(path),
                           'path': path or ''}
                    for name, path in tools.items()}
        info = {
            'pipeline': {
                'name': 'biopytools annorefine (gap-filling)', 'version': mod_ver},
            'phase': 'Phase2 同源查漏补缺|Phase2 homology gap-filling '
                     '(Phase1 BRAKER 由 braker 模块运行|Phase1 BRAKER via braker)',
            'tools': versions,
            'parameters': {
                'gap_min_identity': self.config.gap_min_identity,
                'gap_min_coverage': self.config.gap_min_coverage,
                'gap_min_cds_len': self.config.gap_min_cds_len,
                'require_real_orf': self.config.require_real_orf,
                'gap_coord_zero_overlap': self.config.gap_coord_zero_overlap,
                'unique_reads_only': self.config.unique_reads_only,
                'min_expression_depth': self.config.min_expression_depth,
                'min_coverage_breadth': self.config.min_coverage_breadth,
                'enable_gap_fill': self.config.enable_gap_fill,
                'enable_split': self.config.enable_split,
                'split_min_copy_coverage': self.config.split_min_copy_coverage,
                'threads': self.config.threads,
            },
        }
        out = os.path.join(self.config.pipeline_info_dir, 'software_versions.yml')
        try:
            with open(out, 'w', encoding='utf-8') as f:
                yaml.dump(info, f, default_flow_style=False, allow_unicode=True)
            self.logger.info(f"软件版本信息已保存|Software versions saved: {out}")
        except Exception as e:
            self.logger.warning(f"保存软件版本信息失败|Failed to save versions: {e}")

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
        # enable_small_protein 时 parse 用两通道阈值的较小值, 否则弱命中(id/cov 50/50)在
        # parse 阶段就被 70/80 丢掉, 小蛋白通道拿不到候选; 常规 qc_filter 内部仍按 70/80 严滤
        # |use min(normal,small) at parse so the small lane's weak hits survive; normal
        # qc_filter still applies 70/80 internally, so weak hits can only pass via small lane
        if self.config.enable_small_protein:
            pid = min(self.config.gap_min_identity, self.config.small_min_identity)
            pcov = min(self.config.gap_min_coverage, self.config.small_min_coverage)
        else:
            pid = self.config.gap_min_identity
            pcov = self.config.gap_min_coverage
        hits = parse_miniprot_gff3(self.miniprot_gff, pid, pcov)
        self.logger.info(f"miniprot 命中(过滤后)|hits filtered: {len(hits)}")
        hits = dedupe_hits(hits)   # 全 prot: 同位置多 query 合并|dedup same-locus
        self.logger.info(f"miniprot 命中(去重后)|hits deduped: {len(hits)}")
        return hits

    def run(self):
        """运行查漏补缺流程|Run gap-filling pipeline"""
        self.logger.info("=" * 70)
        self.logger.info("annorefine 查漏补缺: 同源补漏 + 合并拆分 + 质控"
                         "|gap-filling: homology fill + split + QC")
        self.logger.info("=" * 70)

        # Step1 证据扫描|evidence scan
        self._step_evidence_scan()
        if not os.path.exists(self.miniprot_gff):
            self.logger.error("证据扫描未产出|Evidence scan produced nothing")
            return None
        hits = self._load_hits()

        # Step2 漏检/合并分析(路径分治)|gap/merged analysis (path-separated)
        if self.config.skip_gap_analysis:
            self.logger.info("跳过分析|Skipping gap analysis")
            return None
        braker_genes = parse_braker_gff3(self.config.braker_gff3)
        self.logger.info(f"braker 基因|braker genes: {len(braker_genes)}")

        # 2a 纯漏检填补(找全新基因)|gap-fill path (new genes)
        gaps = []
        if self.config.enable_gap_fill:
            gaps = detect_gaps(
                hits, braker_genes, self.config.overlap_cutoff,
                coord_zero_overlap=self.config.gap_coord_zero_overlap)
            self.logger.info(
                f"漏检候选|gap candidates: {len(gaps)} "
                f"(坐标零重叠|coord-zero-overlap={self.config.gap_coord_zero_overlap})")

        # 2b 合并拆分(拆 BRAKER 折叠基因)|merged-split path
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
        self.logger.info(
            f"总候选|total candidates: {len(all_candidate_hits)} "
            f"(gap={len(gaps)}, split={len(merged_hits)})")

        # Step2.5 基因组(ORF检查①) + 表达证据(②④⑤)|genome + expression
        genome = None
        if self.config.require_real_orf:
            genome = load_fasta(self.config.genome)
            self.logger.info(f"基因组加载(ORF检查)|genome loaded: {len(genome)} chroms")
        expression = None
        unique_bam = None
        if self.config.rnaseq_bam:
            unique_bam = prepare_unique_bam(
                self.config.rnaseq_bam, self.config,
                self.cmd_runner, self.logger)
            expression = compute_hit_depth_breadth(
                all_candidate_hits, unique_bam, self.config,
                self.cmd_runner, self.logger)
            if expression is None:
                # depth 计算失败: 醒目告警, 表达过滤将不生效(避免静默丢光所有基因)
                # |depth failed: loud warning, expression filter disabled
                self.logger.error(
                    "表达证据计算失败, 表达过滤②⑤将不生效(仅靠①ORF+③坐标)! "
                    "请检查 BAM/samtools|expression computation FAILED, "
                    "expression filter disabled")
            else:
                self.logger.info(
                    f"表达证据计算|expression computed: {len(expression)} hits")
        else:
            self.logger.warning(
                "无 RNA-seq BAM, 表达过滤②⑤未启用(仅靠①ORF+③坐标)"
                "|no RNA-seq BAM, expression filter off (ORF+coord only)")

        # Step3 质控|QC filter (①ORF + ②⑤表达 + identity/coverage/cds_len/TE)
        repeat_regions = parse_repeat_out(self.config.repeat_out) \
            if self.config.repeat_out else {}
        if not self.config.repeat_out:
            self.logger.warning("未提供 repeat_out, 跳过真TE区排除"
                                "|No repeat_out, skipping real-TE exclusion")

        def _cds_len(h):
            return sum(e - s + 1 for s, e, _ in h.cds_exons)

        # 两候选列表: 常规(>=gap_min_cds_len) + 小蛋白(短且 enable_small_protein)
        # |two candidate lists: normal (>=gap_min_cds_len) + small (short, enable_small_protein)
        normal_candidates = [h for h in all_candidate_hits
                             if _cds_len(h) >= self.config.gap_min_cds_len]
        small_candidates = []
        if self.config.enable_small_protein:
            small_candidates = [h for h in all_candidate_hits
                                if _cds_len(h) < self.config.gap_min_cds_len
                                and _cds_len(h) <= self.config.small_max_cds_len]
            if expression is None:
                self.logger.warning(
                    "小蛋白通道无表达证据 → 退化为 ORF+严同源(70/80)模式|small-protein "
                    "lane has no expression data: degrades to ORF+strict-homology (70/80)")
            self.logger.info(
                f"小蛋白候选|small candidates: {len(small_candidates)} "
                f"(cds_len<{self.config.gap_min_cds_len}bp 且 "
                f"<={self.config.small_max_cds_len}bp)")

        passed = qc_filter(normal_candidates, self.config, repeat_regions,
                           genome=genome, expression=expression)
        self.logger.info(
            f"常规质控通过|normal QC passed: {len(passed)}/{len(normal_candidates)}")
        small_passed = []
        if self.config.enable_small_protein and small_candidates:
            small_passed = qc_filter_small(small_candidates, self.config, repeat_regions,
                                           genome=genome, expression=expression)
            self.logger.info(
                f"小蛋白质控通过|small QC passed: {len(small_passed)}/{len(small_candidates)}")
        passed_ids = {id(h) for h in passed}
        passed_ids.update(id(h) for h in small_passed)

        # ★ 合并拆分门控: 某 merged 基因的 split copies 全未过QC → 不删原 BRAKER 基因(回退保留)
        # |merged-split gating: if no copy passed QC, keep original BRAKER gene
        merged_gene_ids = set()
        for g, copies in merged:
            if any(id(h) in passed_ids for h in copies):
                merged_gene_ids.add(g.gene_id)
        reverted = len(merged) - len(merged_gene_ids)
        if reverted:
            self.logger.info(
                f"合并拆分回退保留|split reverted (copies failed QC): {reverted}")

        # Step3.5 建模型(GFF3): 常规 {prefix}_gap_{N} + 小蛋白 {prefix}_small_gap_{N}
        # |build gene models: normal {prefix}_gap_{N} + small {prefix}_small_gap_{N}
        gap_lines = build_gene_models(passed, self.config.prefix)
        if small_passed:
            small_lines = build_gene_models(small_passed, f"{self.config.prefix}_small")
            gap_lines.extend(small_lines[1:])   # 跳过重复 ##gff-version 头|skip dup header
            self.logger.info(f"小蛋白补基因|small-protein genes added: {len(small_passed)}")
        with open(self.gap_filled_gff3, 'w') as f:
            f.write("\n".join(gap_lines))
            if gap_lines:
                f.write("\n")
        self.logger.info(f"gap_filled 写出|gap_filled written: {self.gap_filled_gff3}")

        # Step4 合并(GFF3)|merge
        if self.config.skip_merge:
            self.logger.info("跳过合并|Skipping merge")
            return self.gap_filled_gff3
        merge_gff3(self.config.braker_gff3, gap_lines,
                   merged_gene_ids, self.merged_gff3)
        self.logger.info(f"merged 写出|merged written: {self.merged_gff3}")

        # Step5 gap 验证报告(复用 expression)|report (reuse expression)
        report_tsv = os.path.join(self.config.gap_dir,
                                  f"{self.config.prefix}.gap_report.tsv")
        from .report import build_gap_report
        build_gap_report(passed, self.config.prefix, self.config.rnaseq_bam,
                         self.config.repeat_out, report_tsv, self.gap_filled_gff3,
                         self.config, self.cmd_runner, self.logger,
                         expression=expression, unique_bam=unique_bam)

        # 清理 unique BAM 临时文件(v2.15: tmp 运行结束清理)|cleanup unique BAM
        if unique_bam:
            for suffix in ('', '.bai', '.csi'):
                p = unique_bam + suffix
                if os.path.exists(p):
                    try:
                        os.remove(p)
                    except OSError:
                        pass

        # 工程收尾: 记录软件版本 + 关键参数(§12.5)|record versions + params
        self._generate_software_versions()

        self.logger.info("=" * 70)
        self.logger.info("annorefine 查漏补缺完成|gap-filling done")
        self.logger.info("=" * 70)
        return self.merged_gff3


def parse_arguments():
    """解析命令行参数(端到端)|Parse CLI arguments (end-to-end)"""
    parser = argparse.ArgumentParser(
        description="annorefine: BRAKER 注释 + 同源查缺补漏 端到端 → 整合 GFF3"
        "|End-to-end: BRAKER + homology gap-filling → integrated GFF3",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例|Example: biopytools annorefine -g genome.fa -s psojae -p prot.fa --rnaseq-dirs r1,r2 -o out/")
    # 必填|Required
    parser.add_argument('-g', '--genome', required=True,
                        help='未mask原始基因组(braker 内部 mask, filling 用未mask)|Unmasked genome')
    parser.add_argument('-s', '--species', required=True,
                        help='物种名(braker 输出命名)|Species name')
    parser.add_argument('-p', '--prot-seq', required=True,
                        help='近缘蛋白(文件或目录, braker+filling 共用)|Protein file/dir')
    parser.add_argument('-o', '--output-dir', required=True,
                        help='输出目录|Output dir')
    # BRAKER 证据|BRAKER evidence
    parser.add_argument('--rnaseq-dirs', help='二代RNA-seq目录(逗号分隔)|RNA-seq dirs')
    parser.add_argument('--isoseq', help='三代转录本(文件或目录)|Iso-seq file/dir')
    # BRAKER 通用|BRAKER general
    parser.add_argument('-t', '--threads', type=int, default=12, help='线程数|Threads (default 12)')
    parser.add_argument('--fungus', action=argparse.BooleanOptionalAction, default=True,
                        help='真菌模式(默认开, --no-fungus 关)|Fungus mode (default on)')
    parser.add_argument('--singularity-image',
                        default='~/software/singularity/braker3_devel.sif',
                        help='Singularity镜像|Singularity image')
    parser.add_argument('--no-singularity', action='store_true',
                        help='不用Singularity|No singularity')
    # BRAKER 步骤|BRAKER steps
    parser.add_argument('--skip-repeat', action='store_true', help='跳过repeat屏蔽|Skip repeat masking')
    parser.add_argument('--skip-repeat-filter', action='store_true',
                        help='跳过repeat库过滤(默认开)|Skip repeat filter')
    parser.add_argument('--skip-rescue', action=argparse.BooleanOptionalAction, default=True,
                        help='跳过证据还原(默认关, --no-skip-rescue 开)|Skip rescue (default on)')
    # 查漏补漏参数|gap-filling params
    parser.add_argument('--split-min-copy-coverage', type=float, default=80,
                        help='保守合并判据:完整拷贝覆盖率%%|Split copy coverage (default 80)')
    parser.add_argument('--no-split', action='store_true', help='关闭合并拆分|Disable merged-gene split')
    parser.add_argument('--repeat-out', help='RepeatMasker .out(默认自动找braker产物)|RepeatMasker out')
    parser.add_argument('--exclude-te-gap', action='store_true',
                        help='质控排除TE区gap(默认不排)|exclude TE-overlap gaps')
    parser.add_argument('--gap-min-identity', type=float, default=70, help='filling identity%%(default 70)')
    parser.add_argument('--gap-min-coverage', type=float, default=80, help='filling coverage%%(default 80)')
    # 通用生物学质控|general bio-QC
    parser.add_argument('--no-real-orf', action='store_true',
                        help='关闭真实完整ORF检查(ATG+stop+3倍数,默认开)|disable real-ORF check (default on)')
    parser.add_argument('--no-coord-zero-overlap', action='store_true',
                        help='关闭gap坐标零重叠(默认开:与BRAKER基因坐标相交不算新基因)|disable coord-zero-overlap (default on)')
    parser.add_argument('--no-unique-reads', action='store_true',
                        help='关闭唯一比对过滤(默认开:多比对reads不算表达)|disable unique-read filter (default on)')
    parser.add_argument('--min-unique-mapq', type=int, default=20,
                        help='唯一比对MAPQ兜底阈值(samtools无-e时)|unique MAPQ fallback (default 20)')
    parser.add_argument('--min-expression-depth', type=float, default=1.0,
                        help='唯一reads平均深度下限(>0)|min unique-read depth (default 1.0)')
    parser.add_argument('--min-coverage-breadth', type=float, default=50.0,
                        help='CDS被唯一reads覆盖广度%%下限|min coverage breadth (default 50)')
    parser.add_argument('--no-gap-fill', action='store_true',
                        help='关闭纯漏检填补(只保留合并拆分)|disable pure gap-fill (split only)')
    # 小蛋白回收通道(通用)|small-protein recovery lane (general)
    parser.add_argument('--recover-small-proteins', action='store_true',
                        help='开启小蛋白回收通道(默认关, 放宽长度找回短蛋白)|enable small-protein lane (default off)')
    parser.add_argument('--small-max-cds-len', type=int, default=450,
                        help='小蛋白CDS上限bp(默认450=150aa)|small max CDS len (default 450)')
    parser.add_argument('--small-min-identity', type=float, default=50.0,
                        help='小蛋白放宽identity%%(默认50, 有表达时)|small min identity (default 50, with expr)')
    parser.add_argument('--small-min-coverage', type=float, default=50.0,
                        help='小蛋白放宽coverage%%(默认50, 有表达时)|small min coverage (default 50, with expr)')
    parser.add_argument('--small-min-expression-depth', type=float, default=1.0,
                        help='小蛋白表达深度下限(默认1.0; effector低表达可调低如0.1)|small min expression depth (default 1.0)')
    parser.add_argument('--small-min-coverage-breadth', type=float, default=60.0,
                        help='小蛋白CDS覆盖广度%%下限(默认60)|small min coverage breadth (default 60)')
    parser.add_argument('--no-small-exclude-te', action='store_true',
                        help='关闭小蛋白TE区排除(默认排; effector常在TE区可关)|disable small-protein TE exclusion (default on)')
    parser.add_argument('--small-strong-homology-identity', type=float, default=95.0,
                        help='强同源直通identity%%阈值(默认95, ≥此值绕过TE/表达过滤)|strong-homology bypass identity (default 95)')
    return parser.parse_args()


def main():
    """主入口|Main entry: BRAKER → 查漏补缺 端到端|BRAKER then gap-filling, end-to-end"""
    args = parse_arguments()

    # 延迟 import(避免 CLI help 时不必要加载)|lazy import
    from ..braker.config import BrakerConfig
    from ..braker.pipeline import BrakerPipeline
    from ..braker.utils import (BrakerLogger, find_protein_files_in_directory,
                                find_long_reads_in_directory)
    from ..braker.main import clean_protein_sequences
    # AnnorefineConfig/Runner 本模块内|local to this module

    # 统一日志(传给 braker+filling, 避免各自重配 root)|unified logger
    log_file = os.path.join(args.output_dir, 'logs', 'annorefine.log')
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    logger = BrakerLogger(log_file).get_logger()

    try:
        logger.info('=' * 70)
        logger.info('annorefine: BRAKER + 查漏补缺 端到端|End-to-end')
        logger.info('=' * 70)

        # 处理 prot_seq(目录识别+清理, braker+filling 共用)|process prot_seq
        prot_seq_file = args.prot_seq
        if prot_seq_file and os.path.isdir(prot_seq_file):
            logger.info('prot_seq 为目录, 自动识别合并|prot_seq is dir, auto-merge')
            prot_seq_file = find_protein_files_in_directory(prot_seq_file, logger=logger)
        if prot_seq_file:
            logger.info('清理蛋白质序列(移除非标准字符)|Cleaning proteins')
            prot_seq_file = clean_protein_sequences(prot_seq_file, logger=logger)

        # 处理 isoseq(目录识别)|process isoseq
        isoseq_file = args.isoseq
        if isoseq_file and os.path.isdir(isoseq_file):
            isoseq_file = find_long_reads_in_directory(isoseq_file, logger=logger)

        # rnaseq_dirs(逗号分隔)|rnaseq_dirs
        rnaseq_dirs = [d.strip() for d in args.rnaseq_dirs.split(',')] if args.rnaseq_dirs else None

        # ===== 阶段1: BRAKER(断点续传, braker.gtf 存在则跳过)|Phase 1: BRAKER =====
        logger.info('-' * 70)
        logger.info('阶段1: BRAKER 注释|Phase 1: BRAKER annotation')
        logger.info('-' * 70)
        bcfg = BrakerConfig(
            genome=args.genome,
            species=args.species,
            prot_seq=prot_seq_file,
            isoseq=isoseq_file,
            rnaseq_dirs=rnaseq_dirs,
            use_singularity=not args.no_singularity,
            singularity_image=args.singularity_image,
            output_dir=args.output_dir,
            threads=args.threads,
            use_fungus=args.fungus,
            skip_repeat=args.skip_repeat,
            skip_repeat_filter=args.skip_repeat_filter,
            skip_rescue=args.skip_rescue,
        )
        bcfg.validate()
        braker_gtf = BrakerPipeline(bcfg, logger).run_pipeline()
        logger.info(f'阶段1 完成, braker.gtf|Phase 1 done: {braker_gtf}')
        # braker 同时输出 braker.gff3(同目录, --gff3), 查漏补缺用 GFF3
        # |braker also outputs braker.gff3 in same dir; gap-filling uses GFF3
        braker_gff3 = braker_gtf.rsplit('.gtf', 1)[0] + '.gff3'
        if not os.path.exists(braker_gff3):
            logger.error(f'braker.gff3 不存在|braker.gff3 not found: {braker_gff3}')
            sys.exit(1)
        logger.info(f'查漏补缺输入 braker.gff3|gap-filling input: {braker_gff3}')

        # ===== 阶段2: 查漏补缺(未mask genome + braker.gff3 + prot)|Phase 2 =====
        # 关键: 用 args.genome(未mask原始), 不是 braker 的 masked genome
        # |uses raw unmasked genome, not braker's masked genome
        logger.info('-' * 70)
        logger.info('阶段2: 查漏补缺|Phase 2: gap-filling')
        logger.info('-' * 70)
        filling_output = os.path.join(args.output_dir, '05_gap_filling')
        # braker 的 RNA-seq BAM(给查漏补缺做表达验证)|braker's RNA-seq BAM for expression
        rnaseq_bam_path = os.path.join(args.output_dir, '03_short_reads',
                                       'rnaseq.sorted.bam')
        rnaseq_bam = [rnaseq_bam_path] if os.path.exists(rnaseq_bam_path) else None
        if rnaseq_bam:
            logger.info(f'表达验证用 RNA-seq BAM|expression BAM: {rnaseq_bam_path}')
        # 自动找 braker 的 RepeatMasker .out(用户 --repeat-out 优先)
        # |auto-find braker's RepeatMasker .out (user --repeat-out takes priority)
        repeat_out = args.repeat_out
        if not repeat_out:
            auto_rep = os.path.join(args.output_dir, '01_repeat_masking',
                                    os.path.basename(args.genome) + '.out')
            if os.path.exists(auto_rep):
                repeat_out = auto_rep
                logger.info(f'自动找到 repeat .out|auto repeat_out: {auto_rep}')
        pcfg = AnnorefineConfig(
            genome=args.genome,
            braker_gff3=braker_gff3,
            prot_seq=prot_seq_file,
            output_dir=filling_output,
            rnaseq_bam=rnaseq_bam,
            threads=args.threads,
            split_min_copy_coverage=args.split_min_copy_coverage,
            enable_split=not args.no_split,
            repeat_out=repeat_out,
            exclude_te_gap=args.exclude_te_gap,
            gap_min_identity=args.gap_min_identity,
            gap_min_coverage=args.gap_min_coverage,
            require_real_orf=not args.no_real_orf,
            gap_coord_zero_overlap=not args.no_coord_zero_overlap,
            unique_reads_only=not args.no_unique_reads,
            min_unique_mapq=args.min_unique_mapq,
            min_expression_depth=args.min_expression_depth,
            min_coverage_breadth=args.min_coverage_breadth,
            enable_gap_fill=not args.no_gap_fill,
            enable_small_protein=args.recover_small_proteins,
            small_max_cds_len=args.small_max_cds_len,
            small_min_identity=args.small_min_identity,
            small_min_coverage=args.small_min_coverage,
            small_min_expression_depth=args.small_min_expression_depth,
            small_min_coverage_breadth=args.small_min_coverage_breadth,
            small_exclude_te=not args.no_small_exclude_te,
            small_strong_homology_identity=args.small_strong_homology_identity,
        )
        result = AnnorefineRunner(pcfg, logger).run()
        logger.info(f'阶段2 完成, 整合 GFF|Phase 2 done: {result}')
        logger.info('=' * 70)
        logger.info('annorefine 端到端完成|annorefine end-to-end done')
        logger.info('=' * 70)
        sys.exit(0)

    except SystemExit:
        raise
    except Exception as e:
        logger.error(f'错误|Error: {e}')
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
