"""
BAM文件批量统计分析命令|BAM File Batch Statistics Analysis Command
"""

import click
import sys
import os


def _lazy_import_bam_stats_main():
    """延迟加载bam_stats主函数|Lazy load bam_stats main function"""
    try:
        from ...bam_stats.main import main as bam_stats_main
        return bam_stats_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_input_exists(input_path):
    """验证输入路径存在|Validate input path exists"""
    if not _is_help_request() and input_path and not os.path.exists(input_path):
        raise click.BadParameter(
            f"输入路径不存在|Input path does not exist: {input_path}"
        )
    return input_path


@click.command(
    short_help='BAM文件批量统计分析|BAM File Batch Statistics',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_input_exists(value) if value else None,
              help='BAM文件或目录|BAM file or directory')
@click.option('--output', '-o',
              default='bam_stats.summary.tsv',
              show_default=True,
              help='输出文件|Output file')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='samtools线程数|Samtools threads')
@click.option('--processes', '-p',
              type=int,
              default=16,
              show_default=True,
              help='并行样本数|Max parallel samples')
@click.option('--reference', '-g',
              help='参考基因组|Reference genome FASTA')
@click.option('--bed-file',
              help='目标区域BED|Target regions BED file')
@click.option('--min-mapq',
              type=int,
              default=20,
              show_default=True,
              help='最小MAPQ|Minimum MAPQ')
@click.option('--max-insert',
              type=int,
              default=1000,
              show_default=True,
              help='最大插入片段|Maximum insert size')
@click.option('--skip-alignment',
              is_flag=True,
              help='跳过比对统计|Skip alignment stats')
@click.option('--skip-coverage',
              is_flag=True,
              help='跳过覆盖度统计|Skip coverage stats')
@click.option('--skip-sequence',
              is_flag=True,
              help='跳过序列特征|Skip sequence features')
@click.option('--skip-insert',
              is_flag=True,
              help='跳过插入片段|Skip insert size stats')
@click.option('--skip-duplicate',
              is_flag=True,
              help='跳过重复统计|Skip duplicate stats')
@click.option('--skip-variation',
              is_flag=True,
              help='跳过变异统计|Skip variation stats')
def bam_stats(input, output, threads, processes,
              reference, bed_file, min_mapq, max_insert,
              skip_alignment, skip_coverage, skip_sequence,
              skip_insert, skip_duplicate, skip_variation):
    """
    BAM文件批量统计分析|BAM File Batch Statistics Analysis

    批量处理BAM文件，输出全局统计长表和染色体级别统计|Batch process BAM files, output summary TSV and per-chromosome TSV

    示例|Examples: biopytools bam-stats -i ./bam_files -o result.summary.tsv
    """
    bam_stats_main = _lazy_import_bam_stats_main()

    args = ['bam_stats.py']
    args.extend(['--input', input])
    args.extend(['--output', output])
    if threads != 24:
        args.extend(['--threads', str(threads)])
    if processes != 16:
        args.extend(['--processes', str(processes)])
    if reference:
        args.extend(['--reference', reference])
    if bed_file:
        args.extend(['--bed-file', bed_file])
    if min_mapq != 20:
        args.extend(['--min-mapq', str(min_mapq)])
    if max_insert != 1000:
        args.extend(['--max-insert', str(max_insert)])
    if skip_alignment:
        args.append('--skip-alignment')
    if skip_coverage:
        args.append('--skip-coverage')
    if skip_sequence:
        args.append('--skip-sequence')
    if skip_insert:
        args.append('--skip-insert')
    if skip_duplicate:
        args.append('--skip-duplicate')
    if skip_variation:
        args.append('--skip-variation')

    original_argv = sys.argv
    sys.argv = args

    try:
        bam_stats_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"运行时错误|Runtime Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
