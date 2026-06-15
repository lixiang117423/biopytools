"""
Wgsim命令|Wgsim Command
"""

import click
import sys
import os


def _lazy_import_wgsim_main():
    """延迟加载wgsim主函数|Lazy load wgsim main function"""
    try:
        from ...wgsim.main import main as wgsim_main
        return wgsim_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在|Validate path exists"""
    if _is_help_request():
        return path
    if not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='Wgsim基因组测序数据模拟|Wgsim genome sequencing simulation',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入基因组文件或目录|Input genome file or directory')
@click.option('-o', '--output-dir',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('-N', '--num-reads',
              type=int,
              default=50000000,
              show_default=True,
              help='模拟reads数量|Number of reads to simulate')
@click.option('-1', '--read-length',
              type=int,
              default=150,
              show_default=True,
              help='reads长度|Read length')
@click.option('-s', '--seed',
              type=int,
              default=0,
              show_default=True,
              help='随机种子|Random seed')
@click.option('-e', '--error-rate',
              type=float,
              default=0.020,
              show_default=True,
              help='测序错误率|Sequencing error rate')
@click.option('-r', '--mutation-rate',
              type=float,
              default=0.001,
              show_default=True,
              help='突变率|Mutation rate')
@click.option('-d', '--outer-distance',
              type=int,
              default=500,
              show_default=True,
              help='外部距离|Outer distance')
@click.option('-D', '--inner-distance',
              type=int,
              default=0,
              show_default=True,
              help='内部距离|Inner distance')
def wgsim(input, output_dir, num_reads, read_length, seed, error_rate, mutation_rate, outer_distance, inner_distance):
    """
    Wgsim基因组测序数据模拟|Wgsim genome sequencing simulation

    示例|Examples: biopytools wgsim -i genome.fna -o output_dir
    """
    wgsim_main = _lazy_import_wgsim_main()

    args = ['wgsim.py']
    args.extend(['-i', input])
    args.extend(['-o', output_dir])

    if num_reads != 50000000:
        args.extend(['-N', str(num_reads)])
    if read_length != 150:
        args.extend(['-1', str(read_length)])
    if seed != 0:
        args.extend(['-s', str(seed)])
    if error_rate != 0.020:
        args.extend(['-e', str(error_rate)])
    if mutation_rate != 0.001:
        args.extend(['-r', str(mutation_rate)])
    if outer_distance != 500:
        args.extend(['-d', str(outer_distance)])
    if inner_distance != 0:
        args.extend(['-D', str(inner_distance)])

    original_argv = sys.argv
    sys.argv = args

    try:
        wgsim_main()
    except SystemExit as e:
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
