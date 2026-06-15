"""
GWAS Lambda GC计算器命令|GWAS Lambda GC Calculator Command
"""

import click
import sys
import os


def _lazy_import_gwas_lambda_main():
    """延迟加载gwas_lambda主函数|Lazy load gwas_lambda main function"""
    try:
        from ...gwas_lambda.main import main as gwas_lambda_main
        return gwas_lambda_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_pattern(pattern):
    """验证搜索模式(仅在非帮助模式)|Validate search pattern (only in non-help mode)"""
    if not _is_help_request():
        if not pattern or not pattern.strip():
            raise click.BadParameter("搜索模式不能为空|Search pattern cannot be empty")
    return pattern.strip() if pattern else pattern


def _validate_threshold(threshold):
    """验证阈值|Validate threshold"""
    if threshold is not None and not (0 < threshold < 1):
        raise click.BadParameter("阈值必须在0和1之间|Threshold must be between 0 and 1")
    return threshold


def _validate_column(column):
    """验证列索引|Validate column index"""
    if column is not None and column < 0:
        raise click.BadParameter("列索引必须为非负数|Column index must be non-negative")
    return column


@click.command(short_help='GWAS Lambda GC计算器|GWAS Lambda GC Calculator',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--pattern', '-p',
              default="feture_*/GWAS_Result.mlm.manht_input",
              show_default=True,
              callback=lambda ctx, param, value: _validate_pattern(value) if value else None,
              help='文件搜索模式|File search pattern')
@click.option('--output', '-o',
              default="Batch_Lambda_Assessment.txt",
              show_default=True,
              help='输出文件名|Output filename')
@click.option('--threshold', '-t',
              type=float,
              default=1e-5,
              show_default=True,
              callback=lambda ctx, param, value: _validate_threshold(value) if value is not None else None,
              help='显著性阈值|Significance threshold')
@click.option('--p-column', '-c',
              type=int,
              default=3,
              show_default=True,
              callback=lambda ctx, param, value: _validate_column(value) if value is not None else None,
              help='P值列索引|P-value column index (0-based)')
@click.option('--output-dir', '-d',
              default="./gwas_lambda_output",
              type=click.Path(),
              show_default=True,
              help='输出目录|Output directory')
def gwas_lambda(pattern, output, threshold, p_column, output_dir):
    """
    GWAS Lambda GC计算器工具|GWAS Lambda GC Calculator Tool

    批量计算GWAS结果的Lambda GC值，用于评估群体分层校正效果
    Calculate Lambda GC values from GWAS results to evaluate population stratification correction

    示例|Examples: biopytools gwas-lambda
    """

    # 延迟加载|Lazy load: import only when actually called
    gwas_lambda_main = _lazy_import_gwas_lambda_main()

    # 构建主函数的参数列表|Build argument list for main function
    args = ['gwas_lambda.py']

    # 添加参数|Add parameters
    args.extend(['--pattern', pattern])
    args.extend(['--output', output])
    args.extend(['--threshold', str(threshold)])
    args.extend(['--p-column', str(p_column)])
    args.extend(['--output-dir', output_dir])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        gwas_lambda_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"执行失败|Execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
