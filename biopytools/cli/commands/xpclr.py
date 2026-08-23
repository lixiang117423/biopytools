"""XP-CLR 跨群体选择信号扫描命令|XP-CLR cross-population selection scan."""
import os
import sys

import click


def _lazy_import_xpclr_main():
    """延迟加载xpclr主函数|Lazy load xpclr main function."""
    try:
        from ...xpclr.main import main as xpclr_main
        return xpclr_main
    except ImportError as e:
        click.echo(f"Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    return any(arg in {'-h', '--help'} for arg in sys.argv)


def _validate_file(path):
    if _is_help_request() or not path:
        return path
    if not os.path.isfile(path):
        raise click.BadParameter(f"文件不存在|File not found: {path}")
    return path


@click.command(
    short_help='XP-CLR跨群体选择信号扫描|XP-CLR cross-population selection scan',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i', required=True,
              callback=lambda ctx, param, value: _validate_file(value),
              help='bgzip VCF(须带.tbi/.csi索引)|bgzipped VCF with tabix index')
@click.option('--samples-a', '-a', required=True,
              callback=lambda ctx, param, value: _validate_file(value),
              help='群体A样本列表(每行一个ID)|Pop A sample list (one ID per line)')
@click.option('--samples-b', '-b', required=True,
              callback=lambda ctx, param, value: _validate_file(value),
              help='群体B样本列表文件|Pop B sample list')
@click.option('--output-dir', '-o', default='xpclr_out', show_default=True,
              help='输出目录|Output directory')
@click.option('--label', default='popA_vs_popB', show_default=True,
              help='结果文件前缀|Result file prefix')
@click.option('--chroms', default=None,
              help='逗号分隔染色体列表(默认VCF全部contig)|Comma-separated chroms (default all)')
@click.option('--size', type=int, default=20000, show_default=True,
              help='窗口大小bp|Window size bp')
@click.option('--step', type=int, default=20000, show_default=True,
              help='滑窗步长bp|Step size bp')
@click.option('--maxsnps', type=int, default=200, show_default=True,
              help='窗口最大SNP数|Max SNPs per window')
@click.option('--minsnps', type=int, default=10, show_default=True,
              help='窗口最小SNP数|Min SNPs per window')
@click.option('--ld', type=float, default=0.95, show_default=True,
              help='LD加权截断|LD cutoff')
@click.option('--phased', is_flag=True, default=False,
              help='数据已phased(更精确r2)|Data phased')
@click.option('--rrate', type=float, default=1e-8, show_default=True,
              help='每碱基重组率|Recombination rate per base')
@click.option('--top-n', type=int, default=50, show_default=True,
              help='Top候选窗口数|Top candidate windows')
@click.option('--backend', type=click.Choice(['xpclrs', 'xpclr']), default='xpclrs',
              show_default=True,
              help='底层工具:xpclrs=Rust高速版(默认),xpclr=python|xpclrs=Rust fast '
                   '(default), xpclr=python')
@click.option('--xpclr-path', default='~/miniforge3/envs/selective_sweep/bin/xpclr',
              show_default=True, help='xpclr(python后端)可执行路径|xpclr (python backend) path')
@click.option('--xpclrs-path', default='~/software/xpclrs/bin/xpclrs',
              show_default=True, help='xpclrs(Rust后端)可执行路径|xpclrs (Rust backend) path')
@click.option('--threads', '-t', type=int, default=12, show_default=True,
              help='线程数(仅xpclrs后端)|Threads, xpclrs backend only')
@click.option('--log-level', default='INFO', show_default=True,
              help='日志级别(DEBUG/INFO/WARNING/ERROR)|Log level')
def xpclr(input, samples_a, samples_b, output_dir, label, chroms, size, step,
          maxsnps, minsnps, ld, phased, rrate, top_n, backend, xpclr_path,
          xpclrs_path, threads, log_level):
    """
    XP-CLR 跨群体选择信号扫描|XP-CLR cross-population selection scan.

    VCF + 群体A/B样本列表 → 逐染色体 XP-CLR → 全基因组合并表 + Top候选窗口表。
    |VCF + pop A/B sample lists -> per-chrom XP-CLR -> genome-wide table + top windows.

    示例|Examples: biopytools xpclr -i pop.vcf.gz -a popA.txt -b popB.txt -o out_dir/
    """

    xpclr_main = _lazy_import_xpclr_main()

    args = ['xpclr.py', '-i', input, '-a', samples_a, '-b', samples_b,
            '-o', output_dir]
    if label != 'popA_vs_popB':
        args.extend(['--label', label])
    if chroms:
        args.extend(['--chroms', chroms])
    if size != 20000:
        args.extend(['--size', str(size)])
    if step != 20000:
        args.extend(['--step', str(step)])
    if maxsnps != 200:
        args.extend(['--maxsnps', str(maxsnps)])
    if minsnps != 10:
        args.extend(['--minsnps', str(minsnps)])
    if ld != 0.95:
        args.extend(['--ld', str(ld)])
    if phased:
        args.append('--phased')
    if rrate != 1e-8:
        args.extend(['--rrate', str(rrate)])
    if top_n != 50:
        args.extend(['--top-n', str(top_n)])
    if backend != 'xpclrs':
        args.extend(['--backend', backend])
    if xpclr_path != '~/miniforge3/envs/selective_sweep/bin/xpclr':
        args.extend(['--xpclr-path', xpclr_path])
    if xpclrs_path != '~/software/xpclrs/bin/xpclrs':
        args.extend(['--xpclrs-path', xpclrs_path])
    if threads != 12:
        args.extend(['--threads', str(threads)])
    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    original_argv = sys.argv
    sys.argv = args
    try:
        xpclr_main()
    except SystemExit as e:
        if e.code not in (0, None):
            raise
    finally:
        sys.argv = original_argv
