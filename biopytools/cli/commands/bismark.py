"""
🧬 Bismark甲基化分析命令 | Bismark Methylation Analysis Command
优化版本：使用懒加载解决响应速度问题
"""

import click
import sys
import os


def _lazy_import_bismark_main():
    """懒加载bismark main函数 | Lazy load bismark main function"""
    try:
        from ...bismark.main import main as bismark_main
        return bismark_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _lazy_import_version():
    """懒加载版本信息 | Lazy load version info"""
    try:
        from ...bismark import __version__
        return __version__
    except ImportError:
        return "1.0.0"  # 默认版本 | Default version


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help', '--version'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path


def _validate_dir_exists(dir_path):
    """验证目录是否存在（仅在非帮助模式下）| Validate directory existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.isdir(dir_path):
        raise click.BadParameter(f"目录不存在 | Directory does not exist: {dir_path}")
    return dir_path


# 加载版本信息用于装饰器
def _get_version_option():
    """获取版本选项装饰器 | Get version option decorator"""
    if _is_help_request():
        # 如果是帮助请求，使用默认版本避免导入
        return click.version_option(version="1.0.0", prog_name='Bismark Pipeline')
    else:
        # 实际运行时才懒加载版本
        version = _lazy_import_version()
        return click.version_option(version=version, prog_name='Bismark Pipeline')


@click.command(short_help="Bismark甲基化分析流程",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
# 版本信息懒加载处理
@click.option('--version', is_flag=True, expose_value=False, is_eager=True,
              callback=lambda ctx, param, value: _handle_version_request() if value else None,
              help='显示版本信息 | Show version information')
# --- Required arguments ---
@click.option('--genome-fa', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='🧬 基因组FASTA文件路径 | Path to genome FASTA file.')
@click.option('--raw-dir', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_dir_exists(value) if value else None,
              help='📁 原始FASTQ数据目录 | Raw FASTQ data directory.')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📂 主输出目录 | Main output directory.')
# --- Optional arguments ---
@click.option('--pattern', '-p',
              default='_1_clean.fq.gz', show_default=True,
              help='🔗 R1文件的后缀模式 | Suffix pattern for R1 files.')
@click.option('--threads', '-j',
              type=int, default=88, show_default=True,
              help='🚀 使用的线程数 | Number of threads to use.')
@click.option('--sort-buffer',
              type=str, default='400G', show_default=True,
              help='💾 提取步骤的排序缓存大小 | Sort buffer size for extraction step.')
@click.option('--include-overlap',
              is_flag=True,
              help='🧪 在提取时【包含】重叠的reads (默认忽略) | Include overlapping reads during extraction.')
def bismark(genome_fa, raw_dir, output_dir, pattern, threads, sort_buffer, include_overlap):
    """🧬 Bismark甲基化分析流程。
    
    一个自动化运行Bismark进行全基因组重亚硫酸盐测序(WGBS)数据分析的流程，
    包括索引构建、比对和甲基化提取。
    
    示例 | Examples:
    
    \b
    # 🎯 基础用法
    biopytools bismark -g genome.fa -r ./raw_data -o ./bismark_results
    
    \b
    # 🚀 使用更多线程和自定义模式
    biopytools bismark -g hg38.fa -r ./cleandata -o ./results -j 96 -p "_1.fq.gz"
        
    \b
    # 🧪 包含重叠的reads
    biopytools bismark -g genome.fa -r ./data -o ./out --include-overlap
    """
    
    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    bismark_main = _lazy_import_bismark_main()
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['bismark.py']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-g', genome_fa])
    args.extend(['-r', raw_dir])
    args.extend(['-o', output_dir])
    
    # 可选参数 ⚙️ | Optional parameters
    if pattern != '_1_clean.fq.gz':
        args.extend(['-p', pattern])
    if threads != 88:
        args.extend(['-j', str(threads)])
    if sort_buffer != '400G':
        args.extend(['--sort-buffer', sort_buffer])
        
    # 处理反向布尔标志 | Handle inverted boolean flag
    # argparse的 '--no-no-overlap' (action='store_false') 意味着默认 no_overlap=True (忽略重叠)
    # 当用户在click中指定 --include-overlap 时，我们必须传递 --no-no-overlap 来让 no_overlap 变为 False
    if include_overlap:
        args.append('--no-no-overlap')

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        bismark_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        if e.code != 0:
            click.secho(f"❌ 脚本执行被终止，退出码: {e.code}", fg='red', err=True)
        sys.exit(e.code)
    except Exception as e:
        click.secho(f"💥 发生未知错误 | An unexpected error occurred: {e}", fg='red', err=True)
        sys.exit(1)
    finally:
        # 无论如何都要恢复原始的 sys.argv | Restore original sys.argv regardless of outcome
        sys.argv = original_argv


def _handle_version_request():
    """处理版本信息请求 | Handle version information request"""
    version = _lazy_import_version()
    click.echo(f"Bismark Pipeline version {version}")
    sys.exit(0)


# 如果直接运行此文件用于测试 | If running this file directly for testing
if __name__ == '__main__':
    bismark()