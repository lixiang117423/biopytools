"""
BAM比对可视化|BAM Alignment Visualization Command
"""

import click
import sys
import os


def _lazy_import_bam_view_main():
    """延迟加载bam_view主函数|Lazy load bam_view main function"""
    try:
        from ...bam_view.main import main as bam_view_main
        return bam_view_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='BAM比对可视化工具|BAM alignment visualization tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--bam', '-b',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='BAM文件路径|BAM file path')
@click.option('--reference', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='参考序列FASTA文件路径|Reference FASTA file path')
@click.option('--region', '-g',
              required=True,
              help='可视化区域(格式: chr:start-end)|Visualization region (format: chr:start-end)')
@click.option('--alignoth-path',
              default='~/miniforge3/envs/alignoth/bin/alignoth',
              show_default=True,
              help='alignoth软件路径|alignoth software path')
@click.option('--output-dir', '-o',
              default='./bam_view_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--output-format', '-f',
              default='html',
              show_default=True,
              type=click.Choice(['html', 'json', 'svg', 'pdf']),
              help='输出格式|Output format')
@click.option('--max-read-depth', '-d',
              type=int,
              default=500,
              show_default=True,
              help='最大reads显示深度|Maximum read depth to display')
@click.option('--max-width', '-w',
              type=int,
              default=1024,
              show_default=True,
              help='最大宽度|Maximum width')
@click.option('--mismatch-display-min-percent',
              type=float,
              default=1.0,
              show_default=True,
              help='显示错配的最小百分比|Minimum percentage of mismatches to display')
@click.option('--vcf', '-v',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF文件路径(高亮变异位点)|VCF file path (highlight variants)')
@click.option('--bed',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='BED文件路径(高亮区域)|BED file path (highlight regions)')
@click.option('--highlight', '-H',
              'highlights',
              multiple=True,
              help='高亮区间(可多次使用, 格式: name:start-end)|Highlight interval (can be used multiple times, format: name:start-end)')
@click.option('--aux-tag', '-x',
              'aux_tags',
              multiple=True,
              help='辅助标签(可多次使用)|Auxiliary tag (can be used multiple times)')
@click.option('--no-embed-js',
              is_flag=True,
              help='不嵌入JavaScript(仅HTML格式)|Do not embed JavaScript (HTML format only)')
@click.option('--plot-all',
              is_flag=True,
              help='绘制所有reads|Plot all reads')
def bam_view(bam, reference, region, alignoth_path, output_dir, output_format,
             max_read_depth, max_width, mismatch_display_min_percent, vcf, bed,
             highlights, aux_tags, no_embed_js, plot_all):
    """
    BAM比对可视化工具|BAM Alignment Visualization Tool

    使用alignoth从BAM文件生成交互式比对可视化|Generate interactive alignment visualization from BAM files using alignoth

    示例|Example: biopytools bam-view -b alignments.bam -r reference.fa -g chr1:1000-2000
    """

    # 延迟加载|Lazy loading
    bam_view_main = _lazy_import_bam_view_main()

    # 构建参数列表|Build argument list
    args = ['bam_view.py']

    # 必需参数|Required parameters
    args.extend(['-b', bam])
    args.extend(['-r', reference])
    args.extend(['-g', region])

    # 可选参数|Optional parameters
    if alignoth_path != '~/miniforge3/envs/alignoth/bin/alignoth':
        args.extend(['--alignoth-path', alignoth_path])

    if output_dir != './bam_view_output':
        args.extend(['-o', output_dir])

    if output_format != 'html':
        args.extend(['-f', output_format])

    if max_read_depth != 500:
        args.extend(['-d', str(max_read_depth)])

    if max_width != 1024:
        args.extend(['-w', str(max_width)])

    if mismatch_display_min_percent != 1.0:
        args.extend(['--mismatch-display-min-percent', str(mismatch_display_min_percent)])

    if vcf:
        args.extend(['-v', vcf])

    if bed:
        args.extend(['--bed', bed])

    for highlight in highlights:
        args.extend(['-H', highlight])

    for aux_tag in aux_tags:
        args.extend(['-x', aux_tag])

    if no_embed_js:
        args.append('--no-embed-js')

    if plot_all:
        args.append('--plot-all')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        bam_view_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
