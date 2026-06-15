"""
HiCanu|HiCanu Genome Assembly Command
"""

import click
import sys
import os


def _lazy_import_hicanu_main():
    """延迟加载hicanu主函数|Lazy load hicanu main function"""
    try:
        from ...hicanu.main import main as hicanu_main
        return hicanu_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
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


@click.command(short_help="HiCanu",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--reads', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入reads文件(FASTA/FASTQ)|Input reads file path (FASTA/FASTQ format)')
@click.option('--genome-size', '-g',
              required=True,
              help='基因组大小(如120m, 1g)|Genome size (e.g., 120m, 1g)')
@click.option('--prefix', '-p',
              required=True,
              help='输出文件前缀|Output file prefix')
@click.option('--output-dir', '-o',
              default='./hicanu_output',
              type=click.Path(),
              show_default=True,
              help='输出目录路径|Output directory path')
@click.option('--canu-path',
              default='~/miniforge3/envs/canu_v.2.3/bin/canu',
              show_default=True,
              help='Canu可执行文件路径|Path to Canu executable')
@click.option('--min-read-length',
              type=int, default=1000,
              show_default=True,
              help='最小read长度|Minimum read length')
@click.option('--min-overlap-length',
              type=int, default=500,
              show_default=True,
              help='最小重叠长度|Minimum overlap length')
@click.option('--corrected-error-rate',
              type=float,
              help='纠错后错误率|Corrected error rate')
@click.option('--raw-error-rate',
              type=float,
              help='原始错误率|Raw error rate')
@click.option('--max-input-coverage',
              type=int,
              help='最大输入覆盖度|Maximum input coverage')
@click.option('--stage',
              type=click.Choice(['haplotype', 'correct', 'trim', 'assemble', 'trim-assemble']),
              default='assemble',
              show_default=True,
              help='组装阶段|Assembly stage')
@click.option('--threads', '-t',
              type=int, default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--memory', '-m',
              default='100G',
              show_default=True,
              help='内存限制|Memory limit')
@click.option('--use-grid',
              is_flag=True,
              help='使用网格引擎|Use grid engine')
@click.option('--grid-options',
              help='网格引擎选项|Grid engine options')
@click.option('--dry-run',
              is_flag=True,
              help='测试运行(不执行)|Dry run (do not execute)')
@click.option('--keep-intermediate',
              is_flag=True,
              help='保留中间文件|Keep intermediate files')
def hicanu(reads, genome_size, prefix, output_dir, canu_path,
           min_read_length, min_overlap_length, corrected_error_rate,
           raw_error_rate, max_input_coverage, stage, threads, memory,
           use_grid, grid_options, dry_run, keep_intermediate):
    """
    HiCanu基因组组装流程|HiCanu Genome Assembly Pipeline

    基于PacBio HiFi reads的基因组组装自动化流程
    Automated genome assembly pipeline based on PacBio HiFi reads

    示例|Examples: biopytools hicanu -i reads.fastq -g 120m -p sample1 -o hicanu_output
    """

    # 延迟加载|Lazy loading: import only when actually called
    hicanu_main = _lazy_import_hicanu_main()

    # 构建主函数参数列表|Build argument list for original main function
    args = ['hicanu.py']

    # 必需参数|Required parameters
    args.extend(['-i', reads])
    args.extend(['-g', genome_size])
    args.extend(['-p', prefix])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if output_dir != './hicanu_output':
        args.extend(['-o', output_dir])

    if canu_path != '~/miniforge3/envs/canu_v.2.3/bin/canu':
        args.extend(['--canu-path', canu_path])

    if min_read_length != 1000:
        args.extend(['--min-read-length', str(min_read_length)])

    if min_overlap_length != 500:
        args.extend(['--min-overlap-length', str(min_overlap_length)])

    if corrected_error_rate is not None:
        args.extend(['--corrected-error-rate', str(corrected_error_rate)])

    if raw_error_rate is not None:
        args.extend(['--raw-error-rate', str(raw_error_rate)])

    if max_input_coverage is not None:
        args.extend(['--max-input-coverage', str(max_input_coverage)])

    if stage != 'assemble':
        args.extend(['--stage', stage])

    if threads != 12:
        args.extend(['-t', str(threads)])

    if memory != '80G':
        args.extend(['-m', memory])

    # 处理选项(布尔标志)|Processing options (boolean flags)
    if use_grid:
        args.append('--use-grid')

    if grid_options:
        args.extend(['--grid-options', grid_options])

    if dry_run:
        args.append('--dry-run')

    if keep_intermediate:
        args.append('--keep-intermediate')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        hicanu_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
