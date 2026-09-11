"""anno_curate命令行入口|anno_curate CLI Entry

示例|Examples: biopytools anno_curate -i anno.gff3 -o out_dir/
"""

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys

from . import __version__
from .config import AnnoCurateConfig
from .diagnosis import (cpc2_diagnosis, diagnosis_fingerprint, merge_diagnosis,
                        rna_coverage_diagnosis, save_fingerprint,
                        should_skip_diagnosis, structural_diagnosis,
                        write_diagnosis)
from .fixer import Gxffixer
from .gff_model import parse_gxf
from .utils import AnnoCurateLogger, format_number


def parse_arguments(argv=None):
    """解析命令行参数|Parse arguments"""
    parser = argparse.ArgumentParser(
        prog='biopytools anno_curate',
        description='GSA提交前注释校正:结构修复+质量诊断(GSAman规则对齐)'
                    '|Pre-GSA annotation curation: structural fix + diagnosis',
        epilog='示例|Examples: biopytools anno_curate -i anno.gff3 -o out_dir/')
    parser.add_argument('-i', '--input', required=True,
                        help='输入GFF3/GTF|Input GFF3/GTF')
    parser.add_argument('-o', '--output-dir', default='./anno_curate_output',
                        help='输出目录|Output directory')
    parser.add_argument('-g', '--genome', default=None,
                        help='基因组FASTA(触发CPC2编码潜能校验)|Genome FASTA '
                             '(enables CPC2 check)')
    parser.add_argument('--rna-bam', nargs='+', default=None,
                        help='转录组BAM(可多个,触发外显子零覆盖校验)|RNA BAMs '
                             '(multiple allowed, enables zero-coverage check)')
    parser.add_argument('--skip-diagnosis', action='store_true', default=False,
                        help='只修不诊|Skip diagnosis')
    parser.add_argument('--no-check-utr', dest='check_utr',
                        action='store_false', default=True,
                        help='关闭UTR比例检查|Disable UTR fraction check')
    parser.add_argument('--min-cds-nt', type=int, default=90,
                        help='CDS过短阈值(nt)|CDS too-short threshold')
    parser.add_argument('--max-cds-nt', type=int, default=10500,
                        help='CDS过长阈值(nt)|CDS too-long threshold')
    parser.add_argument('--utr5-frac-max', type=float, default=0.15,
                        help="5'UTR比例基准|5'UTR fraction baseline")
    parser.add_argument('--utr3-frac-max', type=float, default=0.30,
                        help="3'UTR比例基准|3'UTR fraction baseline")
    parser.add_argument('--relax', type=float, default=0.5,
                        help='松弛系数(有效阈值=基准×(1+relax))|Relax factor')
    parser.add_argument('--cpc2-path', default='',
                        help='CPC2路径(默认查annot域)|CPC2 path '
                             '(defaults to annot domain env)')
    parser.add_argument('--samtools-path', default='',
                        help='samtools路径(默认查align域)|samtools path')
    parser.add_argument('--force', action='store_true', default=False,
                        help='忽略断点全部重跑|Rerun everything')
    parser.add_argument('--log-level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help='日志级别|Log level')
    return parser.parse_args(argv)


def _write_versions(config, logger, out_dir):
    """写 software_versions.yml|Write software_versions.yml"""
    info_dir = os.path.join(out_dir, '00_pipeline_info')
    os.makedirs(info_dir, exist_ok=True)
    lines = [f'anno_curate_module: {__version__}']

    def probe(path, args, logger):
        try:
            # 域环境软件须conda包装(§13);CPC2.py在annot域、samtools在align域
            # |Domain-env tools must be conda-wrapped (§13)
            from biopytools.common.conda_runner import build_conda_command
            cmd = build_conda_command(path, args)
            logger.info(f"命令|Command: {' '.join(cmd)}")
            r = subprocess.run(cmd, shell=False, capture_output=True,
                               text=True, timeout=60)
            text = (r.stdout or '') + (r.stderr or '')
            m = re.search(r'(\d+\.\d+[\w.]*)', text)
            return m.group(1) if m else 'unknown'
        except (subprocess.SubprocessError, OSError):
            return 'unknown'

    if config.genome:
        lines.append(f"cpc2: {probe(config.cpc2_path, ['--version'], logger)}")
    if config.rna_bams:
        lines.append(f"samtools: "
                     f"{probe(config.samtools_path, ['--version'], logger)}")
    path = os.path.join(info_dir, 'software_versions.yml')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')


def run_pipeline(config: AnnoCurateConfig, logger: logging.Logger) -> int:
    """主编排|Main orchestration"""
    sample = config.sample
    fixed_dir = os.path.join(config.output_dir, '01_fixed')
    diag_dir = os.path.join(config.output_dir, '02_diagnosis')
    tmp_dir = os.path.join(config.output_dir, 'tmp')
    os.makedirs(fixed_dir, exist_ok=True)
    # 编排器保证 99_logs 骨架存在(run_pipeline 可被 CLI/测试直接调,
    # 不经 main() 的 AnnoCurateLogger)|Ensure the 99_logs skeleton exists
    os.makedirs(os.path.join(config.output_dir, '99_logs'), exist_ok=True)
    fixed_path = os.path.join(fixed_dir, f'{sample}.fixed.gff3')
    dangling_path = os.path.join(fixed_dir, f'{sample}.fixed.dangling.gff3')
    corrected_path = os.path.join(fixed_dir, f'{sample}_phase_corrected.gff3')
    problematic_path = os.path.join(
        fixed_dir, f'{sample}_phase_problematic.gff3')

    # --force 清掉上一轮 fix 侧旧文件与旧诊断目录:条件附加文件"有内容
    # 才生成",不清理会让本轮未写到的陈旧副本冒充新结果
    # |--force purges last round's fix outputs and stale diagnosis dir;
    # conditional extras are written only when non-empty, so leftovers
    # would masquerade as fresh results
    if config.force:
        for stale in (fixed_path, dangling_path, corrected_path,
                      problematic_path):
            if os.path.exists(stale):
                os.remove(stale)
        shutil.rmtree(diag_dir, ignore_errors=True)

    # ---------- fix 级断点|fix-level resume ----------
    if os.path.exists(fixed_path) and not config.force:
        logger.info(f"修复输出已存在,跳过修复|Fix output exists, skipped: "
                    f"{fixed_path}")
        doc = parse_gxf(fixed_path)
        for v in doc.case_variants:
            logger.warning(v)
    else:
        # 重跑修复(输出缺失):三个条件附加文件以本轮判定为准,先清上一
        # 轮残留,防本轮未写时旧文件冒充结果|Re-fixing: conditional extras
        # reflect THIS round; clear stale copies first so unwritten ones
        # cannot masquerade as fresh output
        for stale in (dangling_path, corrected_path, problematic_path):
            if os.path.exists(stale):
                os.remove(stale)
        doc = parse_gxf(config.input)
        for v in doc.case_variants:
            logger.warning(v)
        logger.info(f"开始修复|Fixing: {format_number(len(doc.records))} "
                    f"records")
        fixer = Gxffixer(doc, logger)
        try:
            result = fixer.run()
        except ValueError as e:
            logger.error(f"修复失败|Fix failed:\n{e}")
            return 1
        written = fixer.write_outputs(fixed_path, dangling_path,
                                      corrected_path, problematic_path)
        logger.info(f"修复完成|Fix done: {result.stats}")
        logger.info(f"附加文件|Extra files: {written}")
        doc = result.doc

    # ---------- 诊断|diagnosis ----------
    if not config.skip_diagnosis:
        tsv = os.path.join(diag_dir, f'{sample}.diagnosis.tsv')
        smy = os.path.join(diag_dir, f'{sample}.diagnosis.summary.txt')
        meta = os.path.join(diag_dir, '.diagnosis_meta.json')
        if should_skip_diagnosis(tsv, meta, config):
            logger.info("诊断输出与参数一致,跳过|Diagnosis unchanged, skipped")
        else:
            os.makedirs(diag_dir, exist_ok=True)
            os.makedirs(tmp_dir, exist_ok=True)
            try:
                issues = structural_diagnosis(doc, config, logger)
                # 裁定A:转录本元数据+全集一遍收集(exon/CDS 不计入 total,
                # cds_len 取真值供 merge_diagnosis 新增 RNA 行用)
                # |Ruling A: collect tx metadata + roster in one pass
                # (exon/CDS never counted in total; real cds_len for
                # merge_diagnosis rows)
                from .gff_model import TRANSCRIPT_RE
                tx_meta = {}
                all_tids = []
                doc.rebuild_indexes()
                for r in doc.records:
                    if TRANSCRIPT_RE.search(r.type):
                        tid = r.attributes.get('ID', '')
                        all_tids.append(tid)
                        cds = [c for c in doc.children_by_parent.get(tid, [])
                               if c.type == 'CDS']
                        cds_len = sum(c.end - c.start + 1 for c in cds)
                        tx_meta[tid] = (r.seqid, r.strand, cds_len)
                extras = {}
                zero_map, no_cov = {}, set()
                if config.genome:
                    nc_tids, extras['cpc2_status'] = cpc2_diagnosis(
                        doc, config, logger, tmp_dir)
                    for tid in nc_tids:
                        issues = _append_tag(issues, tid, 'noncoding', tx_meta)
                else:
                    extras['cpc2_status'] = '跳过|skipped (no --genome)'
                if config.rna_bams:
                    zero_map, no_cov, extras['rna_status'] = \
                        rna_coverage_diagnosis(doc, config, logger, tmp_dir)
                    issues = merge_diagnosis(issues, zero_map, no_cov, tx_meta)
                else:
                    extras['rna_status'] = '跳过|skipped (no --rna-bam)'
                write_diagnosis(issues, all_tids, tsv, smy, extras)
                save_fingerprint(meta, diagnosis_fingerprint(config))
                logger.info(f"诊断完成|Diagnosis done: {len(issues)} issue "
                            f"transcripts")
            finally:
                # 异常也清 tmp(rmtree ignore_errors 已容忍目录不存在)
                # |tmp cleaned even when diagnosis raises (ignore_errors
                # tolerates a missing dir)
                shutil.rmtree(tmp_dir, ignore_errors=True)

    _write_versions(config, logger, config.output_dir)
    return 0


def _append_tag(issues, tid, tag, tx_meta):
    """给已有/新增issue追加标签|Append tag to existing/new issue"""
    from .diagnosis import TranscriptIssue
    for i in issues:
        if i.tid == tid:
            i.tags.append(tag)
            return issues
    seqid, strand, cds_len = tx_meta.get(tid, ('.', '.', 0))
    return issues + [TranscriptIssue(tid, seqid, strand, cds_len, [tag])]


def main():
    """主函数|Main function"""
    args = parse_arguments()
    try:
        config = AnnoCurateConfig(
            input=args.input, output_dir=args.output_dir,
            # argparse dest: --rna-bam → rna_bam(brief 笔误 rna_bams 已修)
            # |dest of --rna-bam is rna_bam (brief typo fixed)
            genome=args.genome, rna_bams=tuple(args.rna_bam or ()),
            skip_diagnosis=args.skip_diagnosis, check_utr=args.check_utr,
            min_cds_nt=args.min_cds_nt, max_cds_nt=args.max_cds_nt,
            utr5_frac_max=args.utr5_frac_max,
            utr3_frac_max=args.utr3_frac_max, relax=args.relax,
            cpc2_path=args.cpc2_path, samtools_path=args.samtools_path,
            force=args.force, log_level=args.log_level)
        config.validate()
        logger = AnnoCurateLogger(
            os.path.join(config.output_dir, '99_logs'),
            config.log_level).get_logger()
        logger.info(f"anno_curate启动|started (input={config.input}, "
                    f"genome={'yes' if config.genome else 'no'}, "
                    f"rna_bams={len(config.rna_bams or ())})")
        sys.exit(run_pipeline(config, logger))
    except ValueError as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
