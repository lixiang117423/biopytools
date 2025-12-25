"""
🧬 GWAS Lambda GC计算命令 | GWAS Lambda GC Calculator Command
优化版本：使用懒加载解决响应速度问题
"""

import click
import sys
import os


def _lazy_import_gwas_lambda_main():
    """懒加载gwas_lambda main函数 | Lazy load gwas_lambda main function"""
    try:
        from ...gwas_lambda.main import main as gwas_lambda_main
        return gwas_lambda_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_pattern(pattern):
    """验证搜索模式（仅在非帮助模式下）| Validate search pattern (only in non-help mode)"""
    if not pattern or not pattern.strip():
        raise click.BadParameter("搜索模式不能为空 | Search pattern cannot be empty")
    return pattern.strip()


def _validate_threshold(threshold):
    """验证阈值 | Validate threshold"""
    if not (0 < threshold < 1):
        raise click.BadParameter("阈值必须在0和1之间 | Threshold must be between 0 and 1")
    return threshold


def _validate_column(column):
    """验证列索引 | Validate column index"""
    if column < 0:
        raise click.BadParameter("列索引必须非负 | Column index must be non-negative")
    return column


@click.command(short_help="GWAS Lambda GC计算工具",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--pattern', '-p',
              default="feture_*/GWAS_Result.mlm.manht_input",
              callback=lambda ctx, param, value: _validate_pattern(value) if value else None,
              help='🔍 文件搜索模式 | File search pattern (default: "feture_*/GWAS_Result.mlm.manht_input")')
@click.option('--output', '-o',
              default="Batch_Lambda_Assessment.txt",
              help='📝 输出文件名 | Output filename (default: "Batch_Lambda_Assessment.txt")')
@click.option('--threshold', '-t',
              type=float,
              default=1e-5,
              callback=lambda ctx, param, value: _validate_threshold(value) if value is not None else None,
              help='🎯 显著性阈值 | Significance threshold (default: 1e-5)')
@click.option('--p-column', '-c',
              type=int,
              default=3,
              callback=lambda ctx, param, value: _validate_column(value) if value is not None else None,
              help='📊 P值所在列索引 | P-value column index (0-based, default: 3)')
@click.option('--output-dir', '-d',
              default="./gwas_lambda_output",
              type=click.Path(),
              help='📁 输出目录 | Output directory (default: "./gwas_lambda_output")')
def gwas_lambda(pattern, output, threshold, p_column, output_dir):
    """
    GWAS Lambda GC计算工具.

    批量分析GWAS结果文件，计算Lambda GC值并评估群体分层情况。
    Lambda GC接近1.0表示结果理想，大于1.1表示可能存在群体分层或假阳性。

    示例 | Examples:

    \b
    # 🎯 使用默认设置分析GWAS结果
    biopytools gwas-lambda

    \b
    # 🔧 指定自定义搜索模式
    biopytools gwas-lambda \\
        --pattern "results/*/gwas.txt" \\
        --output "my_lambda_assessment.txt"

    \b
    # 🎯 调整显著性阈值
    biopytools gwas-lambda \\
        --threshold 1e-6

    \b
    # 📊 指定P值所在列和输出目录
    biopytools gwas-lambda \\
        --p-column 2 \\
        --output-dir "./analysis_output"

    \b
    # 🔍 分析特定格式的文件
    biopytools gwas-lambda \\
        --pattern "chr*/gwas_results.txt" \\
        --output "genome_wide_lambda.txt" \\
        --threshold 5e-8
    """

    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    gwas_lambda_main = _lazy_import_gwas_lambda_main()

    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['gwas_lambda.py']

    # 添加参数 | Add parameters
    if pattern != "feture_*/GWAS_Result.mlm.manht_input":
        args.extend(['--pattern', pattern])
    else:
        args.extend(['--pattern', pattern])

    if output != "Batch_Lambda_Assessment.txt":
        args.extend(['--output', output])
    else:
        args.extend(['--output', output])

    if threshold != 1e-5:
        args.extend(['--threshold', str(threshold)])
    else:
        args.extend(['--threshold', str(threshold)])

    if p_column != 3:
        args.extend(['--p-column', str(p_column)])
    else:
        args.extend(['--p-column', str(p_column)])

    if output_dir != "./gwas_lambda_output":
        args.extend(['--output-dir', output_dir])
    else:
        args.extend(['--output-dir', output_dir])

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 🚀 | Call original main function
        gwas_lambda_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"❌ 错误 | Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv