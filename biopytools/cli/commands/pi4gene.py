"""
基因分组核苷酸多样性计算命令|Nucleotide Diversity per Gene Group Calculation Command
"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...pi4gene.main import main as pi4gene_main
        return pi4gene_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='基因分组核苷酸多样性计算|Calculate nucleotide diversity per gene group',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='序列FASTA文件路径|Input sequence FASTA file path')
@click.option('-d', '--id-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='分组ID文件路径（第一列分组，第二列序列ID）|Group ID file path (col1: group, col2: seq_id)')
@click.option('-o', '--output-dir',
              default='./pi4gene_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--mafft-path',
              default=None,
              help='MAFFT路径|MAFFT path')
def pi4gene(input, id_file, output_dir, threads, mafft_path):
    """
    基因分组核苷酸多样性(pi)计算工具|Nucleotide Diversity (pi) Calculation per Gene Group

    按分组提取序列、MAFFT比对、计算pi
    Extract sequences by group, MAFFT alignment, calculate pi

    示例|Examples: biopytools pi4gene -i genes.fasta -d groups.txt -o pi4gene_output
    """
    pi4gene_main = _lazy_import_main()

    args = ['pi4gene.py']
    args.extend(['-i', input])
    args.extend(['-d', id_file])
    if output_dir != './pi4gene_output':
        args.extend(['-o', output_dir])
    if threads != 12:
        args.extend(['-t', str(threads)])
    if mafft_path:
        args.extend(['--mafft-path', mafft_path])

    original_argv = sys.argv
    sys.argv = args

    try:
        pi4gene_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
