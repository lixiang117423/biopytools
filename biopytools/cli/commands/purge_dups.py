"""
Purge_Dups基因组去冗余|Purge_Dups Genome Deduplication Command
"""

import click
import sys
import os


def _lazy_import_purge_dups_main():
    """延迟加载purge_dups主函数|Lazy load purge_dups main function"""
    try:
        from ...purge_dups.main import main as purge_dups_main
        return purge_dups_main
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
    short_help='Purge_Dups基因组去冗余工具|Purge_Dups genome deduplication tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组组装文件(FASTA格式)|Genome assembly file in FASTA format')
@click.option('--reads', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='测序文件(PacBio/HiFi/illumina)|Sequencing reads file')
@click.option('--output-dir', '-o',
              default='./purge_dups_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--read-type',
              type=click.Choice(['pacbio', 'hifi', 'illumina']),
              default='hifi',
              show_default=True,
              help='测序数据类型|Sequencing data type')
@click.option('--min-fraction',
              type=float,
              default=0.8,
              show_default=True,
              help='最小比例阈值|Minimum fraction threshold')
@click.option('--two-round-chaining',
              is_flag=True,
              default=True,
              help='启用两轮链式匹配|Enable two-round chaining')
@click.option('--no-two-round-chaining',
              is_flag=True,
              help='禁用两轮链式匹配|Disable two-round chaining')
@click.option('--ends-only',
              is_flag=True,
              default=True,
              help='只去除contig末端的冗余|Only remove duplications at contig ends')
@click.option('--no-ends-only',
              is_flag=True,
              help='也去除contig中间的冗余|Also remove duplications in contig middle')
@click.option('--min-primary-length',
              type=int,
              default=10000,
              show_default=True,
              help='最小主contig长度|Minimum primary contig length')
@click.option('--step', '-s',
              type=click.Choice(['1', '2', '3', '4', '5']),
              help='运行指定步骤|Run only specified step:\n'
                   '1: 计算深度|Calculate coverage\n'
                   '2: 计算阈值|Calculate cutoffs\n'
                   '3: 分割比对|Split and align\n'
                   '4: 去冗余|Purge duplications\n'
                   '5: 获取序列|Get sequences')
@click.option('--split-by-n',
              is_flag=True,
              help='split_fa按N分割|split_fa split by N')
@click.option('--manual-cutoffs',
              type=click.Path(exists=True),
              help='手动指定阈值文件|Manual cutoffs file')
def purge_dups(input, reads, output_dir, threads, read_type, min_fraction,
               two_round_chaining, no_two_round_chaining, ends_only, no_ends_only,
               min_primary_length, step, split_by_n, manual_cutoffs):
    """
    Purge_Dups基因组去冗余工具|Purge_Dups Genome Deduplication Tool

    基于测序深度去除基因组组装中的单倍型和重叠序列|
    Remove haplotigs and overlaps in genome assembly based on read depth

    示例|Example: biopytools purge_dups -i assembly.fa -r pacbio_reads.fq -o purged_output
    """

    # 延迟加载|Lazy loading
    purge_dups_main = _lazy_import_purge_dups_main()

    # 构建参数列表|Build argument list
    args = ['purge_dups.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-r', reads])

    # 可选参数|Optional parameters
    if output_dir != './purge_dups_output':
        args.extend(['-o', output_dir])

    if threads != 12:
        args.extend(['-t', str(threads)])

    if read_type != 'hifi':
        args.extend(['--read-type', read_type])

    if min_fraction != 0.8:
        args.extend(['--min-fraction', str(min_fraction)])

    # 处理布尔选项|Handle boolean options
    if two_round_chaining:
        args.append('--two-round-chaining')
    elif no_two_round_chaining:
        args.append('--no-two-round-chaining')

    if ends_only:
        args.append('--ends-only')
    elif no_ends_only:
        args.append('--no-ends-only')

    if min_primary_length != 10000:
        args.extend(['--min-primary-length', str(min_primary_length)])

    if step:
        args.extend(['-s', step])

    if split_by_n:
        args.append('--split-by-n')

    if manual_cutoffs:
        args.extend(['--manual-cutoffs', manual_cutoffs])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        purge_dups_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
