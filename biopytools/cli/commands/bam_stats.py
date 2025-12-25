"""
📊 BAM文件批量统计分析命令 | BAM File Batch Statistics Analysis Command
🚀 高级优化版本：解决--help响应速度问题 ⚡️
"""

import click
import sys
import os


def _lazy_import_bam_stats_main():
    """😴 懒加载bam_stats main函数 | Lazy load bam_stats main function"""
    try:
        from ...bam_stats.main import main as bam_stats_main
        return bam_stats_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """❓ 检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_directory_exists(dir_path):
    """🔍 验证目录是否存在（仅在非帮助模式下）| Validate directory exists (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(dir_path):
        raise click.BadParameter(f"📂❌ 目录不存在 | Directory does not exist: {dir_path}")
    return dir_path


def _validate_parent_dir_writable(file_path):
    """🔍 验证输出文件的父目录是否可写（仅在非帮助模式下）| Validate parent directory is writable (only in non-help mode)"""
    if not _is_help_request():
        parent_dir = os.path.dirname(os.path.abspath(file_path))
        if not os.path.exists(parent_dir):
            try:
                os.makedirs(parent_dir, exist_ok=True)
            except OSError:
                raise click.BadParameter(f"📂❌ 无法创建输出目录 | Cannot create output directory: {parent_dir}")
        elif not os.access(parent_dir, os.W_OK):
            raise click.BadParameter(f"📂❌ 输出目录不可写 | Output directory is not writable: {parent_dir}")
    return file_path


@click.command(short_help="📊 BAM文件批量统计分析",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input-dir', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              help='📁 包含BAM文件的输入目录路径 | Input directory path containing BAM files')
@click.option('--output-file', '-o',
              required=True,
              callback=lambda ctx, param, value: _validate_parent_dir_writable(value) if value else None,
              help='📄 输出报告文件路径 (支持.xlsx和.csv格式) | Output report file path (supports .xlsx and .csv formats)')
@click.option('--processes', '-p',
              type=int,
              help='🔢 并行进程数 | Number of parallel processes (default: all available cores)')
@click.option('--log-dir',
              default='.',
              type=click.Path(),
              help='📁 日志输出目录 | Log output directory (default: current directory)')
def bam_stats(input_dir, output_file, processes, log_dir):
    """
    📊 BAM文件批量统计分析工具.

    💡 功能特性 | Features:

    • 🔄 并行处理多个BAM文件 | Parallel processing of multiple BAM files
    • 📊 提取详细统计信息 | Extract detailed statistics
    • 📈 生成Excel/CSV报告 | Generate Excel/CSV reports
    • 🔍 覆盖度分析 | Coverage analysis
    • 📋 染色体级别统计 | Chromosome-level statistics

    📚 示例 | Examples:

    \b
    # 🎯 基本分析
    biopytools bam-stats -i ./bam_files -o report.xlsx

    \b
    # ⚙️ 指定并行进程数
    biopytools bam-stats -i ./bam_files -o report.csv -p 4

    \b
    # 📁 指定日志目录
    biopytools bam-stats -i ./bam_files -o results/report.xlsx --log-dir ./logs

    \b
    # 📊 生成CSV格式报告
    biopytools bam-stats -i ./bam_files -o analysis_results.csv

    💡 提取的统计信息 | Extracted Statistics:

    • 📊 总reads数 | Total reads
    • 🎯 映射reads数 | Mapped reads
    • 📈 映射率 | Mapping rate
    • 👫 正确配对reads数 | Properly paired reads
    • 📏 覆盖度深度 | Coverage depth
    • 🗺️ 染色体级别统计 | Chromosome-level statistics
    • 📋 覆盖率统计 | Coverage statistics
    """

    # 😴 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    bam_stats_main = _lazy_import_bam_stats_main()

    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['bam_stats.py']

    # 必需参数 📋 | Required parameters
    args.extend(['--input-dir', input_dir])
    args.extend(['--output-file', output_file])

    # 可选参数（只有在非默认值时才添加）⚙️ | Optional parameters (add only when non-default)
    if processes is not None:
        args.extend(['--processes', str(processes)])

    if log_dir != '.':
        args.extend(['--log-dir', log_dir])

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 🚀 | Call original main function
        bam_stats_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"❌ 运行错误 | Runtime Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv