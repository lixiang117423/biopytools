"""
braker4ps 主程序|braker4ps Main Entry
端到端: 先跑 braker 注释(BrakerPipeline), 再接 ps-gene-anno 查漏补缺(PsGeneAnnoRunner)
|End-to-end: run BRAKER then ps-gene-anno gap-filling in one command

约束|Constraint: 不改 braker/ps_gene_anno 源码, 仅 import 调用|import-only
"""

import argparse
import os
import sys


def parse_arguments():
    """解析命令行参数|Parse CLI arguments"""
    parser = argparse.ArgumentParser(
        description="braker4ps: braker 注释 + ps-gene-anno 查漏补缺端到端"
        "|braker + gap-filling end-to-end",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例|Example: biopytools braker4ps -g genome.fa -s psojae -p prot.fa --rnaseq-dirs r1,r2 -o out/")
    # 必填|Required
    parser.add_argument('-g', '--genome', required=True,
                        help='未mask原始基因组(braker 内部 mask, filling 用未mask)|Unmasked genome')
    parser.add_argument('-s', '--species', required=True,
                        help='物种名(braker 输出命名)|Species name')
    parser.add_argument('-p', '--prot-seq', required=True,
                        help='近缘蛋白(文件或目录, braker+filling 共用)|Protein file/dir')
    parser.add_argument('-o', '--output-dir', required=True,
                        help='输出目录(braker 在此, filling 在 output_dir/gap_filling)|Output dir')
    # braker 证据|braker evidence
    parser.add_argument('--rnaseq-dirs', help='二代RNA-seq目录(逗号分隔)|RNA-seq dirs')
    parser.add_argument('--isoseq', help='三代转录本(文件或目录)|Iso-seq file/dir')
    # braker 通用|braker general
    parser.add_argument('-t', '--threads', type=int, default=12, help='线程数|Threads (default 12)')
    parser.add_argument('--fungus', action=argparse.BooleanOptionalAction, default=True,
                        help='真菌模式(疫霉适用, 默认开, --no-fungus 关)|Fungus mode (default on)')
    parser.add_argument('--singularity-image',
                        default='~/software/singularity/braker3_devel.sif',
                        help='Singularity镜像|Singularity image')
    parser.add_argument('--no-singularity', action='store_true',
                        help='不用Singularity|No singularity')
    # braker 步骤|braker steps
    parser.add_argument('--skip-repeat', action='store_true', help='跳过repeat屏蔽|Skip repeat masking')
    parser.add_argument('--skip-repeat-filter', action='store_true',
                        help='跳过repeat库过滤(默认开)|Skip repeat filter')
    parser.add_argument('--skip-rescue', action=argparse.BooleanOptionalAction, default=True,
                        help='跳过证据还原(默认关, --no-skip-rescue 开)|Skip rescue (default on)')
    # filling 参数|filling params
    parser.add_argument('--split-min-copy-coverage', type=float, default=80,
                        help='保守合并判据:完整拷贝覆盖率%|Split copy coverage (default 80)')
    parser.add_argument('--no-split', action='store_true', help='关闭合并拆分|Disable merged-gene split')
    parser.add_argument('--repeat-out', help='RepeatMasker .out(filling真TE排除)|RepeatMasker out')
    parser.add_argument('--exclude-te-gap', action='store_true',
                        help='质控排除TE区gap(默认不排)|exclude TE-overlap gaps')
    parser.add_argument('--gap-min-identity', type=float, default=70, help='filling identity%(default 70)')
    parser.add_argument('--gap-min-coverage', type=float, default=80, help='filling coverage%(default 80)')
    return parser.parse_args()


def main():
    """主入口|Main entry: braker → filling 端到端|braker then filling"""
    args = parse_arguments()

    # 延迟 import(避免 CLI help 时不必要加载)|lazy import
    from ..braker.config import BrakerConfig
    from ..braker.pipeline import BrakerPipeline
    from ..braker.utils import (BrakerLogger, find_protein_files_in_directory,
                                find_long_reads_in_directory)
    from ..braker.main import clean_protein_sequences
    from ..ps_gene_anno.config import PsGeneAnnoConfig
    from ..ps_gene_anno.main import PsGeneAnnoRunner

    # 统一日志(传给 braker+filling, 避免各自重配 root)|unified logger
    log_file = os.path.join(args.output_dir, 'logs', 'braker4ps.log')
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    logger = BrakerLogger(log_file).get_logger()

    try:
        logger.info('=' * 70)
        logger.info('braker4ps: braker + ps-gene-anno 端到端|End-to-end')
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

        # ===== 阶段1: braker(断点续传, braker.gtf 存在则跳过)|Phase 1: braker =====
        logger.info('-' * 70)
        logger.info('阶段1: braker 注释|Phase 1: BRAKER annotation')
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
        # braker 同时输出 braker.gff3(同目录, --gff3), filling 用 GFF3
        # |braker also outputs braker.gff3 in same dir; filling uses GFF3
        braker_gff3 = braker_gtf.rsplit('.gtf', 1)[0] + '.gff3'
        if not os.path.exists(braker_gff3):
            logger.error(f'braker.gff3 不存在|braker.gff3 not found: {braker_gff3}')
            sys.exit(1)
        logger.info(f'filling 输入 braker.gff3|filling input: {braker_gff3}')

        # ===== 阶段2: filling(未mask genome + braker.gtf + prot)|Phase 2: filling =====
        # 关键: filling 用 args.genome(未mask原始), 不是 braker 的 masked genome
        # |filling uses raw unmasked genome, not braker's masked genome
        logger.info('-' * 70)
        logger.info('阶段2: ps-gene-anno 查漏补缺|Phase 2: gap-filling')
        logger.info('-' * 70)
        filling_output = os.path.join(args.output_dir, '05_gap_filling')
        # braker 的 RNA-seq BAM(给 filling gap 报告做表达验证)
        # |braker's RNA-seq BAM for filling gap report expression check
        rnaseq_bam_path = os.path.join(args.output_dir, '03_short_reads',
                                       'rnaseq.sorted.bam')
        rnaseq_bam = [rnaseq_bam_path] if os.path.exists(rnaseq_bam_path) else None
        if rnaseq_bam:
            logger.info(f'gap 报告用 RNA-seq BAM|report BAM: {rnaseq_bam_path}')
        # 自动找 braker 的 RepeatMasker .out(默认, 用户 --repeat-out 优先)
        # |auto-find braker's RepeatMasker .out (user --repeat-out takes priority)
        repeat_out = args.repeat_out
        if not repeat_out:
            auto_rep = os.path.join(args.output_dir, '01_repeat_masking',
                                    os.path.basename(args.genome) + '.out')
            if os.path.exists(auto_rep):
                repeat_out = auto_rep
                logger.info(f'自动找到 repeat .out|auto repeat_out: {auto_rep}')
        pcfg = PsGeneAnnoConfig(
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
        )
        result = PsGeneAnnoRunner(pcfg, logger).run()
        logger.info(f'阶段2 完成, merged.gtf|Phase 2 done: {result}')
        logger.info('=' * 70)
        logger.info('braker4ps 端到端完成|braker4ps end-to-end done')
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
