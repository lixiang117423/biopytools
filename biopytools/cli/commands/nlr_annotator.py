"""
NLR-Annotator分析命令|NLR-Annotator Analysis Command
"""

import click
import sys
import os


def _lazy_import_nlr_annotator_main():
    """延迟加载nlr_annotator主函数|Lazy load nlr_annotator main function"""
    try:
        from ...nlr_annotator.main import main as nlr_annotator_main
        return nlr_annotator_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path exists (only in non-help mode)"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='NLR基因预测工具|NLR gene prediction tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入CDS FASTA文件或目录|Input CDS FASTA file or directory')
@click.option('-o', '--output-dir',
              default='./output',
              help='输出目录|Output directory')
@click.option('-t', '--threads',
              type=int, default=12, show_default=True,
              help='线程数|Number of threads')
@click.option('--sample-suffix',
              default='*.cds.fa', show_default=True,
              help='目录模式下文件匹配后缀|File match suffix for directory mode')
@click.option('--output-gff',
              is_flag=True,
              help='输出GFF文件|Output GFF file')
@click.option('--output-bed',
              is_flag=True,
              help='输出BED文件|Output BED file')
@click.option('--output-motifs',
              is_flag=True,
              help='输出motifs BED文件|Output motifs BED file')
@click.option('--output-alignment',
              is_flag=True,
              help='输出motif比对FASTA|Output motif alignment FASTA')
def nlr_annotator(input, output_dir, threads, sample_suffix, output_gff,
                   output_bed, output_motifs, output_alignment):
    """从CDS序列预测NLR基因|Predict NLR genes from CDS sequences

    示例|Example: biopytools nlr-annotator -i genome.cds.fa -o output_dir/
    """
    nlr_main = _lazy_import_nlr_annotator_main()

    argv = ['nlr_annotator.py']
    argv.extend(['-i', input])
    argv.extend(['-o', output_dir])
    argv.extend(['-t', str(threads)])
    argv.extend(['--sample-suffix', sample_suffix])
    if output_gff:
        argv.append('--output-gff')
    if output_bed:
        argv.append('--output-bed')
    if output_motifs:
        argv.append('--output-motifs')
    if output_alignment:
        argv.append('--output-alignment')

    original_argv = sys.argv
    sys.argv = argv

    try:
        nlr_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
