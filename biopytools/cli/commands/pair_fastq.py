"""
FASTQ配对修复命令|FASTQ Pair Fixing Command
"""

import click
import sys
import os


def _lazy_import_pair_fastq_main():
    """延迟加载pair_fastq主函数|Lazy load pair_fastq main function"""
    try:
        from ...pair_fastq.main import main as pair_fastq_main
        return pair_fastq_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _validate_dir_exists(dir_path):
    """验证目录存在(仅在非帮助模式)|Validate directory exists (only in non-help mode)"""
    help_flags = {'-h', '--help'}
    is_help = any(arg in help_flags for arg in sys.argv)

    if not is_help and not os.path.exists(dir_path):
        raise click.BadParameter(f"目录不存在|Directory does not exist: {dir_path}")
    return dir_path


@click.command(
    short_help='FASTQ配对修复工具|FASTQ pair fixing tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_dir_exists(value) if value else None,
              help='输入目录（包含FASTQ文件）|Input directory containing FASTQ files')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--suffix1',
              default='_1.fq.gz',
              show_default=True,
              help='R1文件后缀|R1 file suffix')
@click.option('--suffix2',
              default='_2.fq.gz',
              show_default=True,
              help='R2文件后缀|R2 file suffix')
@click.option('--tool',
              type=click.Choice(['seqkit', 'repair']),
              default='repair',
              show_default=True,
              help='工具选择|Tool selection (seqkit or repair)')
@click.option('--seqkit-bin',
              default='seqkit',
              show_default=True,
              help='seqkit二进制路径|seqkit binary path')
@click.option('--repair-sh',
              default='repair.sh',
              show_default=True,
              help='repair.sh脚本名称|repair.sh script name')
@click.option('--repair-conda-env',
              default='bbmap_v.39.81',
              show_default=True,
              help='repair.sh的conda环境名称|conda environment for repair.sh')
@click.option('--repair-memory',
              default='300g',
              show_default=True,
              help='repair.sh内存参数|repair.sh memory (-Xmx)')
@click.option('--dry-run',
              is_flag=True,
              help='仅显示命令不执行|Show commands without executing')
@click.option('--verbose',
              is_flag=True,
              help='详细输出|Verbose output')
@click.option('--log-file',
              type=click.Path(),
              help='日志文件路径|Log file path')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
              default='INFO',
              show_default=True,
              help='日志级别|Log level')
def pair_fastq(input, output, threads, suffix1, suffix2, tool, seqkit_bin,
               repair_sh, repair_conda_env, repair_memory, dry_run, verbose, log_file, log_level):
    """
    FASTQ配对修复工具|FASTQ Pair Fixing Tool

    批量修复配对混乱的FASTQ文件|Batch fix paired-end FASTQ files with pairing issues

    示例|Example: biopytools pair-fastq -i raw_data -o fixed_data -t 16
    """

    # 延迟加载|Lazy loading
    pair_fastq_main = _lazy_import_pair_fastq_main()

    # 构建参数列表|Build argument list
    args = ['pair_fastq.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数|Optional parameters
    if threads != 12:
        args.extend(['-t', str(threads)])

    if suffix1 != '_1.fq.gz':
        args.extend(['--suffix1', suffix1])

    if suffix2 != '_2.fq.gz':
        args.extend(['--suffix2', suffix2])

    if tool != 'repair':
        args.extend(['--tool', tool])

    if seqkit_bin != 'seqkit':
        args.extend(['--seqkit-bin', seqkit_bin])

    if repair_sh != 'repair.sh':
        args.extend(['--repair-sh', repair_sh])

    if repair_conda_env != 'bbmap_v.39.81':
        args.extend(['--repair-conda-env', repair_conda_env])

    if repair_memory != '300g':
        args.extend(['--repair-memory', repair_memory])

    # 布尔选项|Boolean options
    if dry_run:
        args.append('--dry-run')

    if verbose:
        args.append('--verbose')

    # 日志选项|Logging options
    if log_file:
        args.extend(['--log-file', log_file])

    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        pair_fastq_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
