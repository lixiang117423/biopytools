"""
WGS混合小种检测命令(单倍体)|WGS mixed-race detection command (haploid).
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
    short_help='WGS混合小种检测(单倍体)|WGS mixed-race detection (haploid)',
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
              help='重复/低复杂度区域BED(可选)|Repeat/low-complexity BED (optional)')
@click.option('--threads', '-t', type=int, default=12, show_default=True, help='线程数|Threads')
@click.option('--kmer-size', '-k', type=int, default=21, show_default=True, help='K-mer大小|K-mer size')
@click.option('--read-length', '-l', type=int, default=150, show_default=True, help='测序读长|Read length')
@click.option('--step', type=int, default=None,
              help='只跑指定步骤1-7(默认全跑)|Run single step 1-7 (default all)')
@click.option('--no-checkpoint', is_flag=True, default=False, help='禁用断点续传|Disable checkpoint')
@click.option('--dry-run', is_flag=True, default=False, help='只打印命令不执行|Print commands only')
@click.option('--min-qual', type=int, default=30, show_default=True, help='变异QUAL下限|Min QUAL')
@click.option('--min-dp', type=int, default=15, show_default=True, help='位点深度下限|Min DP')
@click.option('--min-alt-reads', type=int, default=3, show_default=True, help='ALT支持reads下限|Min ALT reads')
@click.option('--min-coverage', type=int, default=30, show_default=True,
              help='freebayes --min-coverage(默认30)|freebayes min-coverage (default 30)')
@click.option('--min-alt-fraction', type=float, default=0.02, show_default=True,
              help='freebayes --min-alternate-fraction(默认0.02)|freebayes min-alternate-fraction (default 0.02)')
@click.option('--pure-samples', default=None,
              help='已知纯样品(逗号分隔,校准het阈值)|Known-pure samples (comma-sep, calibrate)')
def mixrace(input, clean_fastq_dir, genome, output_dir, repeat_bed, threads, kmer_size,
            read_length, step, no_checkpoint, dry_run, min_qual, min_dp, min_alt_reads,
            min_coverage, min_alt_fraction, pure_samples):
    """
    WGS混合小种检测(单倍体)|WGS mixed-race detection (haploid).

    bwa-mem2+markdup → freebayes -p 1(单倍体)→ 等位频率谱(杂合率/AFS形态/优势株占比)→ 判读。
    |bwa-mem2+markdup -> freebayes -p 1 (haploid) -> AFS (het_rate/shape/dominant) -> verdict.

    示例|Examples: biopytools mixrace -i fastq_dir -g ref.fa -o out_dir/
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
    if min_qual != 30:
        args.extend(['--min-qual', str(min_qual)])
    if min_dp != 15:
        args.extend(['--min-dp', str(min_dp)])
    if min_alt_reads != 3:
        args.extend(['--min-alt-reads', str(min_alt_reads)])
    if min_coverage != 30:
        args.extend(['--min-coverage', str(min_coverage)])
    if min_alt_fraction != 0.02:
        args.extend(['--min-alt-fraction', str(min_alt_fraction)])
    if pure_samples:
        args.extend(['--pure-samples', pure_samples])

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
