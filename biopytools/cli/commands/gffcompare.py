"""
GFFcompare两两比较分析命令|GFFcompare Pairwise Comparison Analysis Command
"""

import click
import sys
import os


def _lazy_import_gffcompare_main():
    """延迟加载gffcompare主函数|Lazy load gffcompare main function"""
    try:
        from ...gffcompare.main import main as gffcompare_main
        return gffcompare_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_input_path(value):
    """验证输入路径存在(仅在非帮助模式)|Validate input path existence (only in non-help mode)"""
    if not _is_help_request() and value and not os.path.exists(value):
        raise click.BadParameter(f"路径不存在|Path does not exist: {value}")
    return value


@click.command(
    short_help='GFF/GTF两两比较分析|GFF/GTF pairwise comparison analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i', multiple=True, required=True,
              callback=_validate_input_path,
              help='输入GFF/GTF文件或文件夹(自动识别)|Input GFF/GTF file(s) or directory (auto-detect)')
@click.option('--output-dir', '-o', default='./gffcompare_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--exon-range', '-e', type=int, default=None,
              help='端部外显子最大允许变异范围|Max terminal exon range')
@click.option('--tss-distance', '-d', type=int, default=None,
              help='转录本起始位点分组距离|TSS grouping distance')
@click.option('-M', '--discard-single-exon-query', is_flag=True,
              help='丢弃单外显子query转录本|Discard single-exon query transcripts')
@click.option('-N', '--discard-single-exon-ref', is_flag=True,
              help='丢弃单外显子reference转录本|Discard single-exon reference transcripts')
@click.option('-R', '--ref-overlap-only', is_flag=True,
              help='仅考虑与query重叠的reference|Only consider reference overlapping query')
@click.option('-Q', '--query-overlap-only', is_flag=True,
              help='仅考虑与reference重叠的query|Only consider query overlapping reference')
@click.option('-T', '--no-tmap-refmap', is_flag=True,
              help='不生成.tmap和.refmap文件|Skip .tmap and .refmap files')
@click.option('--strict-match', is_flag=True,
              help='严格匹配模式|Strict match mode')
@click.option('--cds-match', is_flag=True,
              help='启用CDS链匹配验证|Enable CDS chain matching validation')
@click.option('--genome-seq', '-s', default=None,
              help='基因组序列路径(FASTA)|Genome sequence path (FASTA)')
@click.option('--cprefix', '-p', default=None,
              help='合并GTF中转录本前缀|Transcript prefix in combined GTF')
@click.option('-V', '--verbose-mode', is_flag=True,
              help='详细处理模式|Verbose processing mode')
@click.option('--force', '-f', is_flag=True,
              help='强制重新运行|Force re-run')
@click.option('--gffcompare-path', default=None,
              help='gffcompare软件路径|gffcompare software path')
def gffcompare(input, output_dir, exon_range, tss_distance,
               discard_single_exon_query, discard_single_exon_ref,
               ref_overlap_only, query_overlap_only, no_tmap_refmap,
               strict_match, cds_match, genome_seq, cprefix, verbose_mode,
               force, gffcompare_path):
    """
    GFF/GTF文件两两双向比较分析|GFF/GTF Pairwise Bidirectional Comparison

    使用gffcompare对多个GFF/GTF文件进行两两双向比较，输出分类统计结果
    |Pairwise bidirectional comparison of multiple GFF/GTF files using gffcompare

    示例|Examples: biopytools gffcompare -i sampleA.gff -i sampleB.gtf -o ./output
    示例|Examples: biopytools gffcompare -i ./gff_files/ -o ./output
    """
    gffcompare_main = _lazy_import_gffcompare_main()

    args = ['gffcompare_analysis.py']

    # 必需参数|Required parameters
    args.extend(['-i'] + list(input))

    # 可选参数(仅非默认时添加)|Optional parameters (add only when non-default)
    if output_dir != './gffcompare_output':
        args.extend(['-o', output_dir])

    if exon_range is not None:
        args.extend(['-e', str(exon_range)])

    if tss_distance is not None:
        args.extend(['-d', str(tss_distance)])

    if gffcompare_path:
        args.extend(['--gffcompare-path', gffcompare_path])

    if genome_seq:
        args.extend(['-s', genome_seq])

    if cprefix:
        args.extend(['-p', cprefix])

    # 布尔选项|Boolean options
    if discard_single_exon_query:
        args.append('-M')
    if discard_single_exon_ref:
        args.append('-N')
    if ref_overlap_only:
        args.append('-R')
    if query_overlap_only:
        args.append('-Q')
    if no_tmap_refmap:
        args.append('-T')
    if strict_match:
        args.append('--strict-match')
    if cds_match:
        args.append('--cds-match')
    if verbose_mode:
        args.append('-V')
    if force:
        args.append('-f')

    original_argv = sys.argv
    sys.argv = args

    try:
        gffcompare_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断分析|Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"GFFcompare分析失败|GFFcompare analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
