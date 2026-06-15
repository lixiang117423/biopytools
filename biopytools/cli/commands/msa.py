"""
多序列比对|Multiple Sequence Alignment Command
"""

import click
import sys
import os


def _lazy_import_msa_main():
    """延迟加载MSA主函数|Lazy load MSA main function"""
    try:
        from ...msa.main import main as msa_main
        return msa_main
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
    short_help='多序列比对工具|Multiple sequence alignment tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入序列文件(FASTA格式)|Input sequence file (FASTA format)')
@click.option('--output', '-o',
              required=True,
              help='输出文件前缀|Output file prefix')
@click.option('--method', '-m',
              type=click.Choice(['mafft', 'clustalo', 'muscle', 't_coffee']),
              default='mafft',
              show_default=True,
              help='比对方法|Alignment method')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--format', '-f',
              type=click.Choice(['fasta', 'clustal', 'phylip', 'nexus']),
              default='fasta',
              show_default=True,
              help='输出格式|Output format')
@click.option('--mafft-strategy',
              type=click.Choice(['auto', 'linsi', 'ginsi', 'einsi', 'fftns', 'fftnsi']),
              default='auto',
              show_default=True,
              help='MAFFT比对策略|MAFFT alignment strategy')
@click.option('--mafft-maxiterate',
              type=int,
              default=1000,
              show_default=True,
              help='MAFFT最大迭代次数|MAFFT max iterations')
@click.option('--clustalo-iterations',
              type=int,
              default=0,
              show_default=True,
              help='Clustal Omega迭代次数|Clustal Omega iterations')
@click.option('--muscle-maxiters',
              type=int,
              default=16,
              show_default=True,
              help='MUSCLE最大迭代次数|MUSCLE max iterations')
@click.option('--mafft-path',
              default='mafft',
              show_default=True,
              help='MAFFT程序路径|MAFFT program path')
@click.option('--clustalo-path',
              default='clustalo',
              show_default=True,
              help='Clustal Omega程序路径|Clustal Omega program path')
@click.option('--muscle-path',
              default='muscle',
              show_default=True,
              help='MUSCLE程序路径|MUSCLE program path')
@click.option('--tcoffee-path',
              default='t_coffee',
              show_default=True,
              help='T-Coffee程序路径|T-Coffee program path')
def msa(input, output, method, threads, format,
        mafft_strategy, mafft_maxiterate,
        clustalo_iterations, muscle_maxiters,
        mafft_path, clustalo_path, muscle_path, tcoffee_path):
    """
    多序列比对工具|Multiple Sequence Alignment Tool

    支持MAFFT/Clustal Omega/MUSCLE/T-Coffee进行多序列比对
    Support MAFFT/Clustal Omega/MUSCLE/T-Coffee for multiple sequence alignment

    示例|Examples: biopytools msa -i sequences.fasta -o alignment
    """

    # 延迟加载|Lazy load
    msa_main = _lazy_import_msa_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['msa.py']
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if method != 'mafft':
        args.extend(['-m', method])

    if threads != 88:
        args.extend(['-t', str(threads)])
    if format != 'fasta':
        args.extend(['-f', format])

    # MAFFT参数|MAFFT parameters
    if mafft_strategy != 'auto':
        args.extend(['--mafft-strategy', mafft_strategy])
    if mafft_maxiterate != 1000:
        args.extend(['--mafft-maxiterate', str(mafft_maxiterate)])

    # Clustal Omega参数|Clustal Omega parameters
    if clustalo_iterations != 0:
        args.extend(['--clustalo-iterations', str(clustalo_iterations)])

    # MUSCLE参数|MUSCLE parameters
    if muscle_maxiters != 16:
        args.extend(['--muscle-maxiters', str(muscle_maxiters)])

    # 工具路径|Tool paths
    if mafft_path != 'mafft':
        args.extend(['--mafft-path', mafft_path])
    if clustalo_path != 'clustalo':
        args.extend(['--clustalo-path', clustalo_path])
    if muscle_path != 'muscle':
        args.extend(['--muscle-path', muscle_path])
    if tcoffee_path != 't_coffee':
        args.extend(['--tcoffee-path', tcoffee_path])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        msa_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|MSA interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
