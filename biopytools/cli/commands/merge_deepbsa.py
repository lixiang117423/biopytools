"""
合并Windows版DeepBSA结果|Merge Windows DeepBSA Results
"""

import click
import sys
import os


def _lazy_import_merge_deepbsa_main():
    """延迟加载merge_deepbsa主函数|Lazy load merge_deepbsa main function"""
    try:
        from ...merge_deepbsa.main import main as merge_deepbsa_main
        return merge_deepbsa_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _validate_dir_exists(dir_path):
    """验证目录存在(仅在非帮助模式)|Validate directory existence (only in non-help mode)"""
    help_flags = {'-h', '--help'}
    if not any(arg in help_flags for arg in sys.argv) and dir_path and not os.path.isdir(dir_path):
        raise click.BadParameter(f"目录不存在|Directory does not exist: {dir_path}")
    return dir_path


@click.command(
    short_help='合并Windows版DeepBSA结果|Merge Windows DeepBSA Results',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input-dir',
              required=True,
              callback=lambda ctx, param, value: _validate_dir_exists(value) if value else None,
              help='Windows版DeepBSA输出目录|Windows DeepBSA output directory')
@click.option('-o', '--output-dir',
              required=True,
              help='合并结果输出目录|Merged results output directory')
@click.option('-m', '--methods',
              default=None,
              help='要合并的方法，逗号分隔（默认：全部）|Methods to merge, comma-separated (default: all)')
@click.option('--smooth-func',
              type=click.Choice(['LOWESS', 'Tri-kernel']),
              default='LOWESS',
              show_default=True,
              help='平滑函数（仅fallback时使用）|Smoothing function (only used as fallback)')
@click.option('--smooth-frac',
              type=float,
              default=0.1,
              show_default=True,
              help='平滑窗口比例（仅fallback时使用）|Smoothing window fraction (only used as fallback)')
def merge_deepbsa(input_dir, output_dir, methods, smooth_func, smooth_frac):
    """
    合并Windows版DeepBSA运行结果|Merges Windows DeepBSA results into unified format.

    优先从 npy 文件读取平滑值和原始值（与原始软件完全一致），
    如果 npy 不存在则 fallback 到 values.txt 重新计算。

    示例|Example: biopytools merge-deepbsa -i ./jicaiBSA_Visualize_Results -o ./merged_results
    """

    merge_deepbsa_main = _lazy_import_merge_deepbsa_main()

    args = ['merge_deepbsa.py']
    args.extend(['-i', input_dir])
    args.extend(['-o', output_dir])

    if methods:
        args.extend(['-m', methods])
    if smooth_func != 'LOWESS':
        args.extend(['--smooth-func', smooth_func])
    if smooth_frac != 0.1:
        args.extend(['--smooth-frac', str(smooth_frac)])

    original_argv = sys.argv
    sys.argv = args
    try:
        merge_deepbsa_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
