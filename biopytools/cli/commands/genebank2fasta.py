"""
GenBank序列提取命令|GenBank Sequence Extraction Command
"""

import click
import sys
import os


def _lazy_import_genbank_main():
    """延迟加载GenBank主函数|Lazy load GenBank main function"""
    try:
        from ...genebank2fasta.main import main as genbank_main
        return genbank_main
    except ImportError as e:
        click.echo(f"Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_directory_exists(dir_path):
    """验证目录存在性(仅在非帮助模式下)|Validate directory existence (only in non-help mode)"""
    if not _is_help_request() and dir_path and not os.path.exists(dir_path):
        raise click.BadParameter(f"Directory does not exist: {dir_path}")
    if not _is_help_request() and dir_path and not os.path.isdir(dir_path):
        raise click.BadParameter(f"Path is not a directory: {dir_path}")
    return dir_path


@click.command(
    short_help='GenBank序列提取|GenBank Sequence Extraction',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              help='GenBank文件目录|Input GenBank files directory')
@click.option('--output-dir', '-o',
              default='./genbank_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='并行线程数|Number of parallel threads')
@click.option('--min-length',
              type=int,
              default=10,
              show_default=True,
              help='最小蛋白长度(氨基酸)|Minimum protein length (amino acids)')
@click.option('--phylo', '--create-phylogenetic-matrix',
              is_flag=True,
              help='创建系统发育分析矩阵|Create phylogenetic analysis matrix')
@click.option('--no-sample-sep',
              is_flag=True,
              help='不按样本分离输出|Do not separate output by sample')
@click.option('--no-gene-sep',
              is_flag=True,
              help='不按基因分离输出|Do not separate output by gene')
@click.option('--keep-unknown',
              is_flag=True,
              help='保留未知基因|Keep unknown genes')
def genebank2fasta(input, output_dir, threads, min_length, phylo, no_sample_sep,
                     no_gene_sep, keep_unknown):
    """
    GenBank序列提取工具|GenBank Sequence Extraction Tool

    从GenBank文件提取CDS序列和蛋白序列|Extract CDS and protein sequences from GenBank files

    示例|Examples: biopytools genebank2fasta -i /path/to/genbank -o ./output
    """

    # 延迟加载|Lazy loading: import only when actually called
    genbank_main = _lazy_import_genbank_main()

    # 构建参数列表|Build argument list
    args = ['genebank2fasta.py']

    # Required parameters|必需参数
    args.extend(['-i', input])

    # Optional parameters|可选参数
    if output_dir != './genbank_output':
        args.extend(['-o', output_dir])

    if threads != 88:
        args.extend(['-t', str(threads)])

    if min_length != 10:
        args.extend(['--min-length', str(min_length)])

    # Boolean options|布尔选项
    if phylo:
        args.append('--phylo')

    if no_sample_sep:
        args.append('--no-sample-sep')

    if no_gene_sep:
        args.append('--no-gene-sep')

    if keep_unknown:
        args.append('--keep-unknown')

    # 保存和恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        genbank_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"执行失败|Execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
