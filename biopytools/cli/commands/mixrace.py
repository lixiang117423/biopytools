"""
WGS混合小种检测命令(三分支判读)|WGS mixed-race detection command (three-branch).
"""
import click
import sys
import os


def _lazy_import_mixrace_main():
    """延迟加载mixrace主函数|Lazy load mixrace main function."""
    try:
        from ...mixrace.main import main as mixrace_main
        return mixrace_main
    except ImportError as e:
        click.echo(f"Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    return any(arg in {'-h', '--help'} for arg in sys.argv)


def _validate_dir(path):
    if _is_help_request() or not path:
        return path
    if not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


def _validate_file(path):
    if _is_help_request() or not path:
        return path
    if not os.path.isfile(path):
        raise click.BadParameter(f"文件不存在|File not found: {path}")
    return path


@click.command(
    short_help='WGS混合小种检测(三分支判读)|WGS mixed-race detection (three-branch)',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i', default=None,
              callback=lambda ctx, param, value: _validate_dir(value) if value else None,
              help='原始FASTQ目录(与--clean-fastq-dir二选一)|Raw FASTQ dir (or --clean-fastq-dir)')
@click.option('--clean-fastq-dir', default=None,
              callback=lambda ctx, param, value: _validate_dir(value) if value else None,
              help='已清洗FASTQ目录(给则跳过QC)|Clean FASTQ dir (skip QC)')
@click.option('--genome', '-g', required=True,
              callback=lambda ctx, param, value: _validate_file(value) if value else None,
              help='参考基因组FASTA|Reference genome FASTA')
@click.option('--output-dir', '-o', default='mixrace_out', show_default=True,
              help='输出目录|Output directory')
@click.option('--repeat-bed', default=None,
              callback=lambda ctx, param, value: _validate_file(value) if value else None,
              help='额外排除区域BED(与自动热点并集)|Extra exclude BED (merged with hotspots)')
@click.option('--host-genome', default=None,
              callback=lambda ctx, param, value: _validate_file(value) if value else None,
              help='寄主基因组FASTA(给则比对寄主并整对剔除寄主reads,报告寄主占比)'
              '|Host genome FASTA (deplete host reads, report host rate)')
@click.option('--min-mapq', type=int, default=20, show_default=True,
              help='比对质量阈值:mapped reads提取+统计口径(0=不过滤)'
              '|Min MAPQ: mapped-read extraction + stats (0=off)')
@click.option('--threads', '-t', type=int, default=12, show_default=True, help='线程数|Threads')
@click.option('--kmer-size', '-k', type=int, default=21, show_default=True, help='K-mer大小|K-mer size')
@click.option('--read-length', '-l', type=int, default=150, show_default=True, help='测序读长|Read length')
@click.option('--step', type=int, default=None,
              help='只跑指定步骤1-5(1=QC+寄主剔除 2=GTX 3=评估判读 4=k-mer 5=图+报告)'
              '|Run single step 1-5 (default all)')
@click.option('--no-checkpoint', is_flag=True, default=False, help='禁用断点续传|Disable checkpoint')
@click.option('--dry-run', is_flag=True, default=False, help='只打印命令不执行|Print commands only')
@click.option('--pure-het-threshold', type=float, default=0.001, show_default=True,
              help='总杂合率低于此值判纯菌(0.001=0.1%)|Pure threshold')
@click.option('--partner-alt-rate', type=float, default=0.8, show_default=True,
              help='混合伴侣:ALT携带率阈值|Partner ALT-carrier threshold')
@click.option('--partner-hom-rate', type=float, default=0.5, show_default=True,
              help='混合伴侣:伴侣纯合1/1占比阈值|Partner homozygous threshold')
@click.option('--min-sites', type=int, default=1000, show_default=True,
              help='最低有GT位点数,低于判uncertain|Min called sites')
@click.option('--window-size', type=int, default=100000, show_default=True,
              help='热点窗口大小bp|Hotspot window size')
@click.option('--hotspot-fold', type=float, default=2.0, show_default=True,
              help='热点:窗口杂合率>该倍数x自身全基因组率|Hotspot fold')
@click.option('--hotspot-min-median', type=float, default=0.10, show_default=True,
              help='热点:窗口候选中位杂合率下限|Hotspot min median rate')
def mixrace(input, clean_fastq_dir, genome, output_dir, repeat_bed, host_genome, min_mapq,
            threads, kmer_size, read_length, step, no_checkpoint, dry_run,
            pure_het_threshold, partner_alt_rate, partner_hom_rate, min_sites,
            window_size, hotspot_fold, hotspot_min_median):
    """
    WGS混合小种检测(三分支判读)|WGS mixed-race detection (three-branch).

    fastp(+可选寄主剔除) → GTX联合calling → 四层杂合评估(L1 AD/DP排错 →
    L2 shared/private+混合伴侣 → L3 窗口 → L4 热点排除)+群体结构 → 三分支判读
    (纯菌/优势菌株/混杂菌株)+实验建议 → mapped reads k-mer → 全套图+报告。
    |QC(+host depletion) -> GTX joint calling -> four-layer het evaluation +
    partner analysis -> three-branch verdict + advice -> mapped-read k-mer -> report.

    示例|Examples: biopytools mixrace -i fastq_dir -g ref.fa --host-genome host.fa -o out_dir/
    """

    mixrace_main = _lazy_import_mixrace_main()

    args = ['mixrace.py']
    if clean_fastq_dir:
        args.extend(['--clean-fastq-dir', clean_fastq_dir])
    elif input:
        args.extend(['-i', input])
    args.extend(['-g', genome, '-o', output_dir])
    if repeat_bed:
        args.extend(['--repeat-bed', repeat_bed])
    if host_genome:
        args.extend(['--host-genome', host_genome])
    if min_mapq != 20:
        args.extend(['--min-mapq', str(min_mapq)])
    if threads != 12:
        args.extend(['-t', str(threads)])
    if kmer_size != 21:
        args.extend(['-k', str(kmer_size)])
    if read_length != 150:
        args.extend(['-l', str(read_length)])
    if step is not None:
        args.extend(['--step', str(step)])
    if no_checkpoint:
        args.append('--no-checkpoint')
    if dry_run:
        args.append('--dry-run')
    if pure_het_threshold != 0.001:
        args.extend(['--pure-het-threshold', str(pure_het_threshold)])
    if partner_alt_rate != 0.8:
        args.extend(['--partner-alt-rate', str(partner_alt_rate)])
    if partner_hom_rate != 0.5:
        args.extend(['--partner-hom-rate', str(partner_hom_rate)])
    if min_sites != 1000:
        args.extend(['--min-sites', str(min_sites)])
    if window_size != 100000:
        args.extend(['--window-size', str(window_size)])
    if hotspot_fold != 2.0:
        args.extend(['--hotspot-fold', str(hotspot_fold)])
    if hotspot_min_median != 0.10:
        args.extend(['--hotspot-min-median', str(hotspot_min_median)])

    original_argv = sys.argv
    sys.argv = args
    try:
        mixrace_main()
    except SystemExit as e:
        if e.code not in (0, None):
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"执行失败|Execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
