"""aliner命令|aliner Command (a-liner synteny visualization pipeline)"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载aliner主函数|Lazy load aliner main"""
    try:
        from ...aliner.main import main as aliner_main
        return aliner_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否帮助请求|Check help request"""
    return any(a in ('-h', '--help') for a in sys.argv)


def _validate_path(path):
    """验证路径存在|Validate path exists"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help="a-liner共线性可视化(FASTA->minimap2->图)|a-liner synteny (FASTA->minimap2->plot)",
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--ref', 'ref', required=True,
              callback=lambda c, p, v: _validate_path(v),
              help='参考基因组FASTA|Reference genome FASTA')
@click.option('--query', 'query', required=True,
              callback=lambda c, p, v: _validate_path(v),
              help='查询基因组FASTA|Query genome FASTA')
@click.option('--ref-seqs', 'ref_seqs', required=True,
              help='ref侧序列(逗号分隔,如 chrZ,chrW 或 chrZ:1-30000000)|ref-side seqs')
@click.option('--query-seqs', 'query_seqs', required=True,
              help='query侧序列(逗号分隔,与ref等长)|query-side seqs')
@click.option('-o', '--output-dir', 'output_dir', default='./aliner_output',
              help='输出目录|Output directory')
@click.option('--out-prefix', 'out_prefix', default='synteny', help='输出前缀|Output prefix')
@click.option('--preset', type=click.Choice(['asm5', 'asm10', 'asm20']), default='asm5',
              help='minimap2预设(近缘asm5)|minimap2 preset')
@click.option('--min-identity', 'min_identity', type=int, default=70, help='identity阈值%%|identity threshold')
@click.option('--min-alignment-len', 'min_alignment_len', type=int, default=1000,
              help='最小比对长度|min alignment length')
@click.option('--colormap', type=click.Choice(['0', '1', '2', '3', '4', '5']), default='5', help='配色|colormap')
@click.option('--figure-size', 'figure_size', type=float, nargs=2, default=[6, 0],
              metavar='W H', help='图尺寸[宽 高]|figure size')
@click.option('-t', '--threads', 'threads', type=int, default=12, help='线程数|Threads')
@click.option('--extra-args', 'extra_args', default='', help='透传a-liner参数|pass-through args')
def aliner(ref, query, ref_seqs, query_seqs, output_dir, out_prefix, preset,
           min_identity, min_alignment_len, colormap, figure_size, threads, extra_args):
    """FASTA->minimap2->a-liner 共线性可视化pipeline"""

    aliner_main = _lazy_import_main()

    args = ['aliner.py',
            '--ref', ref, '--query', query,
            '--ref-seqs', ref_seqs, '--query-seqs', query_seqs,
            '-o', output_dir, '--out-prefix', out_prefix,
            '--preset', preset,
            '--min-identity', str(min_identity),
            '--min-alignment-len', str(min_alignment_len),
            '--colormap', colormap,
            '--figure-size', str(figure_size[0]), str(figure_size[1]),
            '--threads', str(threads)]
    if extra_args:
        args.extend(['--extra-args', extra_args])

    original_argv = sys.argv
    sys.argv = args
    try:
        aliner_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
