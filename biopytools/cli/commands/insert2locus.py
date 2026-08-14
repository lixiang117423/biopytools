"""insert2locus CLI包装器|insert2locus CLI wrapper"""

import os
import sys

import click


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...insert2locus.main import main as module_main
        return module_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在|Validate path exists"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(short_help="转基因插入位点提取(步移+完整locus+验证)"
                          "|Transgenic insertion locus extraction")
@click.option('-i', '--input', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='fastq目录或R1文件(自动识别样本)|fastq dir or R1 file (auto-detect samples)')
@click.option('-f', '--insert-fasta', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='插入序列fasta(载体+片段)|Insert fasta (vector+fragment)')
@click.option('-o', '--output-dir', required=True,
              help='输出目录|Output directory')
@click.option('-t', '--threads', default=12, show_default=True,
              help='线程数|Threads')
@click.option('--sort-mem', default='2G', show_default=True,
              help='samtools sort单线程内存|samtools sort per-thread memory')
@click.option('--read1-suffix', default=None,
              help='R1后缀(默认自动检测)|R1 suffix (auto-detect by default)')
@click.option('--max-rounds', default=30, show_default=True,
              type=int, help='步移最大轮数|Max walking rounds')
@click.option('--min-softclip', default=25, show_default=True, type=int,
              help='诱饵softclip最短长度|Min softclip for bait')
@click.option('--min-unmapped', default=400, show_default=True, type=int,
              help='未比对contig作诱饵最短长度|Min unmapped length for bait')
@click.option('--min-growth', default=50, show_default=True, type=int,
              help='收敛判定增量阈值|Growth threshold for convergence')
@click.option('--mapq-min', default=1, show_default=True, type=int,
              help='招募最低MAPQ|Min MAPQ for recruitment')
@click.option('--repeat-cap', default=10000, show_default=True, type=int,
              help='单轮新招reads上限(撞重复区)|Per-round new-read cap')
@click.option('--junction-flank', default=50, show_default=True, type=int,
              help='边界报告最短侧翼|Min flank for border report')
@click.option('--tdna-fasta', default=None,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else value,
              help='单独插入序列fasta(区分insert与载体骨架)|Standalone T-DNA fasta')
@click.option('--target-flank', default=2000, show_default=True, type=int,
              help='LB/RB目标侧翼长度|Target LB/RB flank length')
@click.option('--force', is_flag=True, default=False,
              help='忽略断点全部重跑|Ignore checkpoints and rerun')
@click.option('--log-level', default='INFO', show_default=True,
              help='日志级别|Log level')
def insert2locus(input, insert_fasta, output_dir, threads, sort_mem,
                 read1_suffix, max_rounds, min_softclip, min_unmapped,
                 min_growth, mapq_min, repeat_cap, junction_flank, tdna_fasta,
                 target_flank, force, log_level):
    """
    转基因插入位点提取:soft-clip钓取+迭代步移+完整locus重构+WGS覆盖验证+HTML报告|
    Transgenic insertion locus extraction: soft-clip fishing + iterative walking
    + complete locus reconstruction + coverage verification + HTML report

    示例|Examples: biopytools insert2locus -i fq_dir/ -f insert.fasta -o output/
    """

    module_main = _lazy_import_main()

    args = ['insert2locus.py']
    args.extend(['-i', input])
    args.extend(['-f', insert_fasta])
    args.extend(['-o', output_dir])
    if threads != 12:
        args.extend(['-t', str(threads)])
    if sort_mem != '2G':
        args.extend(['--sort-mem', sort_mem])
    if read1_suffix is not None:
        args.extend(['--read1-suffix', read1_suffix])
    if max_rounds != 30:
        args.extend(['--max-rounds', str(max_rounds)])
    if min_softclip != 25:
        args.extend(['--min-softclip', str(min_softclip)])
    if min_unmapped != 400:
        args.extend(['--min-unmapped', str(min_unmapped)])
    if min_growth != 50:
        args.extend(['--min-growth', str(min_growth)])
    if mapq_min != 1:
        args.extend(['--mapq-min', str(mapq_min)])
    if repeat_cap != 10000:
        args.extend(['--repeat-cap', str(repeat_cap)])
    if junction_flank != 50:
        args.extend(['--junction-flank', str(junction_flank)])
    if tdna_fasta is not None:
        args.extend(['--tdna-fasta', tdna_fasta])
    if target_flank != 2000:
        args.extend(['--target-flank', str(target_flank)])
    if force:
        args.append('--force')
    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    original_argv = sys.argv
    sys.argv = args

    try:
        module_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
