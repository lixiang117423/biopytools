"""
WGS混合小种检测命令|WGS mixed-race detection command.
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
    """检查是否为帮助请求|Check if this is a help request."""
    return any(arg in {'-h', '--help'} for arg in sys.argv)


def _validate_dir(input_path):
    """验证输入目录存在|Validate input directory exists."""
    if _is_help_request() or not input_path:
        return input_path
    if not os.path.exists(input_path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {input_path}")
    return input_path


def _validate_file(path):
    """验证文件存在|Validate file exists."""
    if _is_help_request() or not path:
        return path
    if not os.path.isfile(path):
        raise click.BadParameter(f"文件不存在|File not found: {path}")
    return path


@click.command(
    short_help='WGS混合小种检测|WGS mixed-race detection',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_dir(value) if value else None,
              help='输入FASTQ目录(自动配对)|Input FASTQ dir (auto-pair)')
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file(value) if value else None,
              help='参考基因组FASTA|Reference genome FASTA')
@click.option('--output-dir', '-o',
              default='mixrace_out',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--repeat-bed',
              default=None,
              callback=lambda ctx, param, value: _validate_file(value) if value else None,
              help='重复/低复杂度区域BED(可选,跳过去重)|Repeat/low-complexity BED (optional)')
@click.option('--threads', '-t', type=int, default=12, show_default=True,
              help='线程数|Number of threads')
@click.option('--kmer-size', '-k', type=int, default=21, show_default=True,
              help='K-mer大小|K-mer size')
@click.option('--read-length', '-l', type=int, default=150, show_default=True,
              help='测序读长|Read length')
@click.option('--step', type=int, default=None,
              help='只跑指定步骤1-8(默认全跑)|Run single step 1-8 (default all)')
@click.option('--no-checkpoint', is_flag=True, default=False,
              help='禁用断点续传|Disable checkpoint resume')
@click.option('--dry-run', is_flag=True, default=False,
              help='只打印命令不执行|Print commands without executing')
@click.option('--min-qual', type=int, default=30, show_default=True,
              help='变异QUAL下限|Min variant QUAL')
@click.option('--min-dp', type=int, default=15, show_default=True,
              help='位点深度下限|Min site DP')
@click.option('--min-alt-reads', type=int, default=3, show_default=True,
              help='ALT等位支持reads下限|Min ALT-supporting reads')
@click.option('--vaf-mid-low', type=float, default=0.05, show_default=True,
              help='VAF中间频率占比-纯阈值|VAF mid-ratio pure threshold')
@click.option('--vaf-mid-high', type=float, default=0.15, show_default=True,
              help='VAF中间频率占比-杂阈值|VAF mid-ratio mixed threshold')
@click.option('--multiallelic-low', type=float, default=0.005, show_default=True,
              help='多等位占比-纯阈值|Multiallelic ratio pure threshold')
@click.option('--multiallelic-high', type=float, default=0.03, show_default=True,
              help='多等位占比-杂阈值|Multiallelic ratio mixed threshold')
@click.option('--fws-cutoff', type=float, default=0.95, show_default=True,
              help='Fws下限(单基因型)|Fws cutoff (single genotype)')
@click.option('--min-depth', type=float, default=15.0, show_default=True,
              help='判读可信平均深度下限|Min mean depth for confident verdict')
@click.option('--pure-samples', default=None,
              help='已知纯样品(逗号分隔,用于校准阈值)|Known-pure samples (comma-sep, calibrate)')
def mixrace(input, genome, output_dir, repeat_bed, threads, kmer_size, read_length,
            step, no_checkpoint, dry_run, min_qual, min_dp, min_alt_reads,
            vaf_mid_low, vaf_mid_high, multiallelic_low, multiallelic_high,
            fws_cutoff, min_depth, pure_samples):
    """
    WGS混合小种检测|WGS mixed-race detection.

    输入fastq目录+参考基因组,输出每样品"疑似纯/疑似混合/灰色"判读+依据。
    |Input fastq dir + reference genome; outputs per-sample pure/mixed/uncertain verdict + rationale.

    示例|Examples: biopytools mixrace -i fastq_dir -g ref.fa -o out_dir/
    """

    mixrace_main = _lazy_import_mixrace_main()

    args = ['mixrace.py']
    # 必需参数|required
    args.extend(['-i', input, '-g', genome, '-o', output_dir])
    # 可选参数(非默认才传)|optional (only when non-default)
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
    if vaf_mid_low != 0.05:
        args.extend(['--vaf-mid-low', str(vaf_mid_low)])
    if vaf_mid_high != 0.15:
        args.extend(['--vaf-mid-high', str(vaf_mid_high)])
    if multiallelic_low != 0.005:
        args.extend(['--multiallelic-low', str(multiallelic_low)])
    if multiallelic_high != 0.03:
        args.extend(['--multiallelic-high', str(multiallelic_high)])
    if fws_cutoff != 0.95:
        args.extend(['--fws-cutoff', str(fws_cutoff)])
    if min_depth != 15.0:
        args.extend(['--min-depth', str(min_depth)])
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
