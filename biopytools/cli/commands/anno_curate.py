"""anno_curate 注释校正命令|anno_curate annotation curation command"""

import os
import sys

import click


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...anno_curate.main import main as anno_curate_main
        return anno_curate_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    return any(arg in {'-h', '--help'} for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(非帮助模式)|Validate path exists (non-help mode)"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='GSA提交前注释校正(修复+诊断)|Pre-GSA annotation curation '
               '(fix + diagnosis)',
    context_settings=dict(help_option_names=['-h', '--help'],
                          max_content_width=120))
@click.option('--input', '-i', default=None, required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              type=click.Path(), help='输入GFF3/GTF|Input GFF3/GTF')
@click.option('--output-dir', '-o', default='./anno_curate_output',
              show_default=True, type=click.Path(),
              help='输出目录|Output directory')
@click.option('--genome', '-g', default=None,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              type=click.Path(),
              help='基因组FASTA(触发CPC2)|Genome FASTA (enables CPC2)')
@click.option('--rna-bam', multiple=True,
              callback=lambda ctx, param, value: [
                  _validate_path_exists(v) for v in value],
              help='转录组BAM,可多次传(触发零覆盖校验)|RNA BAM, '
                   'repeatable (enables zero-coverage check)')
@click.option('--skip-diagnosis', is_flag=True, default=False,
              show_default=True, help='只修不诊|Skip diagnosis')
@click.option('--check-utr/--no-check-utr', default=True, show_default=True,
              help='UTR比例检查|UTR fraction check')
@click.option('--min-cds-nt', default=90, show_default=True, type=int)
@click.option('--max-cds-nt', default=10500, show_default=True, type=int)
@click.option('--utr5-frac-max', default=0.15, show_default=True, type=float)
@click.option('--utr3-frac-max', default=0.30, show_default=True, type=float)
@click.option('--relax', default=0.5, show_default=True, type=float,
              help='松弛系数,有效阈值=基准×(1+relax)|Relax factor')
@click.option('--cpc2-path', default=None, help='CPC2路径|CPC2 path')
@click.option('--samtools-path', default=None, help='samtools路径|samtools')
@click.option('--force', is_flag=True, default=False,
              help='忽略断点重跑|Rerun all')
@click.option('--log-level', default='INFO', show_default=True,
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']))
def anno_curate(input, output_dir, genome, rna_bam, skip_diagnosis,
                check_utr, min_cds_nt, max_cds_nt, utr5_frac_max,
                utr3_frac_max, relax, cpc2_path, samtools_path, force,
                log_level):
    """GSA提交前注释校正:结构修复+质量诊断(GSAman规则对齐)
    |Pre-GSA annotation curation: structural fix + diagnosis

    示例|Examples: biopytools anno_curate -i anno.gff3 -o out_dir/
    """
    anno_curate_main = _lazy_import_main()
    args = ['anno_curate', '-i', input, '-o', output_dir]
    if genome:
        args.extend(['-g', genome])
    if rna_bam:
        args.extend(['--rna-bam'] + list(rna_bam))
    if skip_diagnosis:
        args.append('--skip-diagnosis')
    if not check_utr:
        args.append('--no-check-utr')
    args.extend(['--min-cds-nt', str(min_cds_nt),
                 '--max-cds-nt', str(max_cds_nt),
                 '--utr5-frac-max', str(utr5_frac_max),
                 '--utr3-frac-max', str(utr3_frac_max),
                 '--relax', str(relax)])
    if cpc2_path:
        args.extend(['--cpc2-path', cpc2_path])
    if samtools_path:
        args.extend(['--samtools-path', samtools_path])
    if force:
        args.append('--force')
    args.extend(['--log-level', log_level])
    original_argv = sys.argv
    sys.argv = args
    try:
        anno_curate_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
