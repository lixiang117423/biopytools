"""
annorefine: BRAKER + 查漏补缺 端到端 → 整合 GFF3|braker + gap-filling end-to-end
"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载主函数|Lazy load main"""
    try:
        from ...annorefine import main as anno_module
        return anno_module.main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    return any(a in {'-h', '--help'} for a in sys.argv)


def _validate_file(path):
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"文件不存在|File not found: {path}")
    return path


@click.command(
    short_help='BRAKER+查漏补缺端到端→整合GFF3|braker + gap-filling end-to-end',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-g', '--genome', required=True,
              callback=lambda c, p, v: _validate_file(v) if v else None,
              help='未mask原始基因组|Unmasked genome')
@click.option('-s', '--species', required=True, help='物种名|Species name')
@click.option('-p', '--prot-seq', required=True,
              callback=lambda c, p, v: _validate_file(v) if v else None,
              help='近缘蛋白(文件或目录)|Protein file/dir')
@click.option('-o', '--output-dir', required=True, help='输出目录|Output dir')
@click.option('--rnaseq-dirs', help='二代RNA-seq目录(逗号分隔)|RNA-seq dirs')
@click.option('--isoseq', help='三代转录本(文件或目录)|Iso-seq file/dir')
@click.option('-t', '--threads', type=int, default=12, show_default=True, help='线程数|Threads')
@click.option('--fungus/--no-fungus', default=True, show_default=True,
              help='真菌模式(疫霉适用)|Fungus mode')
@click.option('--no-singularity', is_flag=True, help='不用Singularity|No singularity')
@click.option('--skip-repeat', is_flag=True, help='跳过repeat屏蔽|Skip repeat masking')
@click.option('--skip-repeat-filter', is_flag=True, help='跳过repeat库过滤|Skip repeat filter')
@click.option('--skip-rescue/--no-skip-rescue', default=True, show_default=True,
              help='证据还原(默认关)|Rescue (default off)')
@click.option('--split-min-copy-coverage', type=float, default=80, show_default=True,
              help='保守合并判据:完整拷贝覆盖率%|Split copy coverage')
@click.option('--no-split', is_flag=True, help='关闭合并拆分|Disable split')
@click.option('--repeat-out', help='RepeatMasker .out(filling真TE排除)|RepeatMasker out')
@click.option('--exclude-te-gap', is_flag=True, help='质控排除TE区gap(默认不排)|exclude TE-overlap gaps')
@click.option('--no-real-orf', is_flag=True,
              help='关闭真实完整ORF检查(默认开)|disable real-ORF check (default on)')
@click.option('--no-coord-zero-overlap', is_flag=True,
              help='关闭gap坐标零重叠(默认开)|disable coord-zero-overlap (default on)')
@click.option('--no-unique-reads', is_flag=True,
              help='关闭唯一比对过滤(默认开)|disable unique-read filter (default on)')
@click.option('--min-unique-mapq', type=int, default=20, show_default=True,
              help='唯一比对MAPQ兜底阈值|unique MAPQ fallback')
@click.option('--min-expression-depth', type=float, default=1.0, show_default=True,
              help='唯一reads平均深度下限|min unique-read depth')
@click.option('--min-coverage-breadth', type=float, default=50.0, show_default=True,
              help='CDS覆盖广度%下限|min coverage breadth')
@click.option('--no-gap-fill', is_flag=True,
              help='关闭纯漏检填补(只保留合并拆分)|disable pure gap-fill (split only)')
@click.option('--recover-small-proteins', is_flag=True,
              help='开启小蛋白回收通道(默认关, 通用)|enable small-protein lane (default off)')
@click.option('--small-max-cds-len', type=int, default=450, show_default=True,
              help='小蛋白CDS上限bp|small max CDS len')
@click.option('--small-min-identity', type=float, default=50.0, show_default=True,
              help='小蛋白放宽identity%(有表达时)|small min identity')
@click.option('--small-min-coverage', type=float, default=50.0, show_default=True,
              help='小蛋白放宽coverage%(有表达时)|small min coverage')
@click.option('--small-min-expression-depth', type=float, default=1.0, show_default=True,
              help='小蛋白表达深度下限(effector低表达可调低如0.1)|small min expression depth')
@click.option('--small-min-coverage-breadth', type=float, default=60.0, show_default=True,
              help='小蛋白CDS覆盖广度%下限|small min coverage breadth')
@click.option('--no-small-exclude-te', is_flag=True,
              help='关闭小蛋白TE区排除(effector常在TE区可关)|disable small-protein TE exclusion')
@click.option('--small-strong-homology-identity', type=float, default=95.0, show_default=True,
              help='强同源直通identity%阈值(≥此值绕过TE/表达过滤)|strong-homology bypass identity')
def annorefine(genome, species, prot_seq, output_dir, rnaseq_dirs, isoseq,
              threads, fungus, no_singularity, skip_repeat, skip_repeat_filter,
              skip_rescue, split_min_copy_coverage, no_split, repeat_out,
              exclude_te_gap, no_real_orf, no_coord_zero_overlap, no_unique_reads,
              min_unique_mapq, min_expression_depth, min_coverage_breadth, no_gap_fill,
              recover_small_proteins, small_max_cds_len, small_min_identity, small_min_coverage,
              small_min_expression_depth, small_min_coverage_breadth, no_small_exclude_te,
              small_strong_homology_identity):
    """
    BRAKER 注释 + 查漏补缺 端到端 → 整合 GFF3|braker + gap-filling end-to-end

    示例|Example: biopytools annorefine -g genome.fa -s psojae -p prot.fa -o out/
    """
    anno_main = _lazy_import_main()
    args = ['annorefine.py', '-g', genome, '-s', species, '-p', prot_seq, '-o', output_dir]
    if rnaseq_dirs:
        args.extend(['--rnaseq-dirs', rnaseq_dirs])
    if isoseq:
        args.extend(['--isoseq', isoseq])
    if threads != 12:
        args.extend(['-t', str(threads)])
    if not fungus:
        args.append('--no-fungus')
    if no_singularity:
        args.append('--no-singularity')
    if skip_repeat:
        args.append('--skip-repeat')
    if skip_repeat_filter:
        args.append('--skip-repeat-filter')
    if not skip_rescue:
        args.append('--no-skip-rescue')
    if split_min_copy_coverage != 80:
        args.extend(['--split-min-copy-coverage', str(split_min_copy_coverage)])
    if no_split:
        args.append('--no-split')
    if repeat_out:
        args.extend(['--repeat-out', repeat_out])
    if exclude_te_gap:
        args.append('--exclude-te-gap')
    if no_real_orf:
        args.append('--no-real-orf')
    if no_coord_zero_overlap:
        args.append('--no-coord-zero-overlap')
    if no_unique_reads:
        args.append('--no-unique-reads')
    if min_unique_mapq != 20:
        args.extend(['--min-unique-mapq', str(min_unique_mapq)])
    if min_expression_depth != 1.0:
        args.extend(['--min-expression-depth', str(min_expression_depth)])
    if min_coverage_breadth != 50.0:
        args.extend(['--min-coverage-breadth', str(min_coverage_breadth)])
    if no_gap_fill:
        args.append('--no-gap-fill')
    if recover_small_proteins:
        args.append('--recover-small-proteins')
    if small_max_cds_len != 450:
        args.extend(['--small-max-cds-len', str(small_max_cds_len)])
    if small_min_identity != 50.0:
        args.extend(['--small-min-identity', str(small_min_identity)])
    if small_min_coverage != 50.0:
        args.extend(['--small-min-coverage', str(small_min_coverage)])
    if small_min_expression_depth != 1.0:
        args.extend(['--small-min-expression-depth', str(small_min_expression_depth)])
    if small_min_coverage_breadth != 60.0:
        args.extend(['--small-min-coverage-breadth', str(small_min_coverage_breadth)])
    if no_small_exclude_te:
        args.append('--no-small-exclude-te')
    if small_strong_homology_identity != 95.0:
        args.extend(['--small-strong-homology-identity', str(small_strong_homology_identity)])

    original_argv = sys.argv
    sys.argv = args
    try:
        anno_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
