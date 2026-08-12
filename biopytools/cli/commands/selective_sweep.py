"""
选择性扫荡检测命令|Selective Sweep Detection Command
"""

import click
import sys
import os


def _lazy_import_selective_sweep_main():
    """延迟加载selective_sweep主函数|Lazy load selective_sweep main function"""
    try:
        from ...selective_sweep.main import main as selective_sweep_main
        return selective_sweep_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在|Validate file existence"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(
            f"文件不存在|File does not exist: {file_path}"
        )
    return file_path


@click.command(
    short_help='选择性扫荡检测|Selective sweep detection',
    context_settings=dict(
        help_option_names=['-h', '--help'], max_content_width=120
    )
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: (
                  _validate_file_exists(value) if value else None
              ),
              help='输入VCF文件|Input VCF file')
@click.option('--pop-info', '-p',
              required=True,
              callback=lambda ctx, param, value: (
                  _validate_file_exists(value) if value else None
              ),
              help='群体分组文件(样品ID<TAB>分组,无表头)|Population info file')
@click.option('--output-dir', '-o',
              default='./selective_sweep_output',
              type=click.Path(),
              show_default=True,
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--win',
              type=int,
              default=50000,
              show_default=True,
              help='窗口大小|Window size')
@click.option('--step',
              type=int,
              default=50000,
              show_default=True,
              help='窗口步长|Window step')
@click.option('--top-quantile',
              type=float,
              default=0.01,
              show_default=True,
              help='候选阈值分位数|Candidate threshold quantile')
@click.option('--merge-gap',
              type=int,
              default=100000,
              show_default=True,
              help='候选窗口合并最大间隔|Max gap for merging candidate windows')
@click.option('--min-maf',
              type=float,
              default=0.05,
              show_default=True,
              help='过滤MAF阈值|Filter MAF threshold')
@click.option('--max-missing',
              type=float,
              default=0.10,
              show_default=True,
              help='过滤缺失率阈值|Filter missing rate threshold')
@click.option('--raisd-window',
              type=int,
              default=50,
              show_default=True,
              help='RAiSD SNP窗口|RAiSD SNP window')
@click.option('--raisd-min-samples',
              type=int,
              default=15,
              show_default=True,
              help='低样本量阈值(低于此值MU分量默认排除)|Low sample threshold (MU excluded below)')
@click.option('--include-mu-low-n',
              is_flag=True,
              default=False,
              help='低样本群体也加入MU分量|Include MU component for low-n pops')
@click.option('--sweed-grid',
              type=int,
              default=10000,
              show_default=True,
              help='SweeD CLR网格点数|SweeD CLR grid points')
@click.option('--sweed-min-samples',
              type=int,
              default=15,
              show_default=True,
              help='低样本量阈值(低于此值CLR分量默认排除)|Low sample threshold (CLR excluded below)')
@click.option('--include-sweed-low-n',
              is_flag=True,
              default=False,
              help='低样本群体也加入CLR分量|Include CLR component for low-n pops')
@click.option('--sweed-unfolded',
              is_flag=True,
              default=False,
              help='使用未折叠SFS(需祖先状态,默认折叠)|Unfolded SFS (requires ancestral state; folded by default)')
@click.option('--xpclr-maxsnps',
              type=int,
              default=200,
              show_default=True,
              help='XP-CLR窗口最大SNP数|XP-CLR max SNPs per window')
@click.option('--xpclr-minsnps',
              type=int,
              default=10,
              show_default=True,
              help='XP-CLR窗口最小SNP数|XP-CLR min SNPs per window')
@click.option('--xpclr-ld',
              type=float,
              default=0.95,
              show_default=True,
              help='XP-CLR LD加权截断|XP-CLR LD cutoff')
@click.option('--xpclr-min-samples',
              type=int,
              default=15,
              show_default=True,
              help='低样本量阈值(任一群体低于此值XP-CLR分量默认排除)|Low sample threshold (XP-CLR excluded if either pop below)')
@click.option('--include-xpclr-low-n',
              is_flag=True,
              default=False,
              help='低样本群体对也加入XP-CLR分量|Include XP-CLR component for low-n pairs')
@click.option('--log-level',
              default='INFO',
              show_default=True,
              help='日志级别|Log level')
def selective_sweep(input, pop_info, output_dir, threads, win, step,
                    top_quantile, merge_gap, min_maf, max_missing,
                    raisd_window, raisd_min_samples, include_mu_low_n,
                    sweed_grid, sweed_min_samples, include_sweed_low_n,
                    sweed_unfolded, xpclr_maxsnps, xpclr_minsnps, xpclr_ld,
                    xpclr_min_samples, include_xpclr_low_n, log_level):
    """
    选择性扫荡检测工具|Selective Sweep Detection Tool

    输入VCF与群体分组信息,自动完成过滤/拆分/统计计算/合并打分,输出候选扫荡区域
    |Input VCF + pop info; auto filter/split/calc/merge-score; output candidate regions

    示例|Examples: biopytools selective-sweep -i variants.vcf.gz -p pops.txt -o sweep_output
    """

    selective_sweep_main = _lazy_import_selective_sweep_main()

    args = ['selective_sweep.py']
    args.extend(['-i', input])
    args.extend(['-p', pop_info])

    if output_dir != './selective_sweep_output':
        args.extend(['-o', output_dir])
    if threads != 12:
        args.extend(['-t', str(threads)])
    if win != 50000:
        args.extend(['--win', str(win)])
    if step != 50000:
        args.extend(['--step', str(step)])
    if top_quantile != 0.01:
        args.extend(['--top-quantile', str(top_quantile)])
    if merge_gap != 100000:
        args.extend(['--merge-gap', str(merge_gap)])
    if min_maf != 0.05:
        args.extend(['--min-maf', str(min_maf)])
    if max_missing != 0.10:
        args.extend(['--max-missing', str(max_missing)])
    if raisd_window != 50:
        args.extend(['--raisd-window', str(raisd_window)])
    if raisd_min_samples != 15:
        args.extend(['--raisd-min-samples', str(raisd_min_samples)])
    if include_mu_low_n:
        args.extend(['--include-mu-low-n'])
    if sweed_grid != 10000:
        args.extend(['--sweed-grid', str(sweed_grid)])
    if sweed_min_samples != 15:
        args.extend(['--sweed-min-samples', str(sweed_min_samples)])
    if include_sweed_low_n:
        args.extend(['--include-sweed-low-n'])
    if sweed_unfolded:
        args.extend(['--sweed-unfolded'])
    if xpclr_maxsnps != 200:
        args.extend(['--xpclr-maxsnps', str(xpclr_maxsnps)])
    if xpclr_minsnps != 10:
        args.extend(['--xpclr-minsnps', str(xpclr_minsnps)])
    if xpclr_ld != 0.95:
        args.extend(['--xpclr-ld', str(xpclr_ld)])
    if xpclr_min_samples != 15:
        args.extend(['--xpclr-min-samples', str(xpclr_min_samples)])
    if include_xpclr_low_n:
        args.extend(['--include-xpclr-low-n'])
    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    original_argv = sys.argv
    sys.argv = args

    try:
        selective_sweep_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
