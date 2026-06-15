"""
系统发育树构建|Phylogenetic Tree Builder Command
"""

import click
import sys
import os


def _lazy_import_phylo_main():
    """延迟加载系统发育树构建主函数|Lazy load phylo main function"""
    try:
        from ...mafft_fasttree.main import main as phylo_main
        return phylo_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
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
    short_help='系统发育树构建工具|Phylogenetic tree builder tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入序列文件(FASTA格式)|Input sequence file (FASTA format)')
@click.option('--output', '-o',
              default='./phylo_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--seq-type',
              type=click.Choice(['protein', 'nucleotide']),
              help='序列类型(未指定时自动检测)|Sequence type (auto-detect if not specified)')
@click.option('--mafft-params',
              default='--auto',
              show_default=True,
              help='MAFFT额外参数|Additional MAFFT parameters')
@click.option('--fasttree-params',
              default='',
              show_default=True,
              help='FastTree额外参数|Additional FastTree parameters')
@click.option('--mafft-path',
              default='mafft',
              show_default=True,
              help='MAFFT软件路径|MAFFT software path')
@click.option('--fasttree-path',
              default='fasttree',
              show_default=True,
              help='FastTree软件路径|FastTree software path')
def mafft_fasttree(input, output, threads, seq_type, mafft_params, fasttree_params,
                   mafft_path, fasttree_path):
    """
    系统发育树构建工具|Phylogenetic Tree Builder Tool

    使用MAFFT进行多序列比对,FastTree构建系统发育树
    Multiple sequence alignment with MAFFT, phylogenetic tree construction with FastTree

    示例|Examples: biopytools mafft-fasttree -i sequences.fa -o phylo_results
    """

    # 延迟加载|Lazy load
    phylo_main = _lazy_import_phylo_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['mafft_fasttree.py']
    args.extend(['-i', input])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if output != './phylo_output':
        args.extend(['-o', output])
    if threads != 88:
        args.extend(['-t', str(threads)])
    if seq_type:
        args.extend(['--seq-type', seq_type])
    if mafft_params != '--auto':
        args.extend(['--mafft-params', mafft_params])
    if fasttree_params:
        args.extend(['--fasttree-params', fasttree_params])
    if mafft_path != 'mafft':
        args.extend(['--mafft-path', mafft_path])
    if fasttree_path != 'fasttree':
        args.extend(['--fasttree-path', fasttree_path])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        phylo_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Pipeline interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
