"""FASTA序列长度统计 CLI 包装|seq-len CLI wrapper"""

import os
import sys

import click


def _lazy_import_seq_len_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...seq_len.main import main as seq_len_main
        return seq_len_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    return any(arg in {'-h', '--help'} for arg in sys.argv)


def _validate_path_exists(path):
    """验证文件或文件夹存在(仅在非帮助模式)|Validate file/folder exists (non-help mode only)"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='FASTA序列长度统计(文件/文件夹,合并+汇总)|FASTA length stats (file/folder, merged + summary)',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i', required=True,
              callback=lambda c, p, v: _validate_path_exists(v),
              help='FASTA 文件或文件夹|FASTA file or folder')
@click.option('--output', '-o', required=True,
              help='输出 TSV 路径或目录|Output TSV path or directory')
@click.option('--prefix', default=None,
              help='输出前缀(默认取输入名)|Output prefix (default: input name)')
@click.option('--min-length', type=int, default=0, show_default=True,
              help='最小长度过滤(0=不过滤)|Min length filter (0=no filter)')
@click.option('--sort', 'sort', is_flag=True,
              help='按长度降序(默认保持输入顺序)|Sort by length descending')
@click.option('--no-summary', 'no_summary', is_flag=True,
              help='不输出汇总表|Skip summary table')
@click.option('--log-file', default=None, help='日志文件|Log file')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
              default='INFO', show_default=True, help='日志级别|Log level')
@click.option('--verbose', '-v', is_flag=True, help='详细日志|Verbose')
def seq_len(input, output, prefix, min_length, sort, no_summary, log_file, log_level, verbose):
    """
    FASTA序列长度统计|FASTA sequence length statistics

    输入单个 FASTA 文件或一个含 FASTA 的文件夹,计算每条序列长度并输出 TSV;
    文件夹模式合并所有文件并加 source_file 列,同时输出 N50 等汇总表。
    |Input a FASTA file or a folder of FASTAs; compute per-sequence length to TSV.
    Folder mode merges all files with a source_file column and writes an N50 summary.

    示例|Examples: biopytools seq-len -i genome.fa -o out.tsv
    """

    seq_len_main = _lazy_import_seq_len_main()

    args = ['seq_len.py', '-i', input, '-o', output]
    if prefix:
        args.extend(['--prefix', prefix])
    if min_length != 0:
        args.extend(['--min-length', str(min_length)])
    if sort:
        args.append('--sort')
    if no_summary:
        args.append('--no-summary')
    if log_file:
        args.extend(['--log-file', log_file])
    if log_level != 'INFO':
        args.extend(['--log-level', log_level])
    if verbose:
        args.append('--verbose')

    original_argv = sys.argv
    sys.argv = args
    try:
        seq_len_main()
    except SystemExit as e:
        if e.code not in (0, None):
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
