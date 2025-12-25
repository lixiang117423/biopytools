"""
🧬 甲烷循环基因丰度分析命令 | Methane Cycle Gene Abundance Analysis Command
优化版本：使用懒加载解决响应速度问题
"""

import click
import sys
import os


def _lazy_import_mcyc_main():
    """懒加载mcyc main函数 | Lazy load mcyc main function"""
    try:
        from ...mcyc.main import main as mcyc_main
        return mcyc_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path


def _validate_threads(threads):
    """验证线程数 | Validate thread count"""
    if threads < 1:
        raise click.BadParameter("线程数必须大于0 | Thread count must be greater than 0")
    return threads


@click.command(short_help="甲烷循环基因丰度分析工具",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.argument('input-list',
                callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
                required=True)
@click.option('--output', '-o',
              default=None,
              type=click.Path(),
              help='📁 输出目录路径 | Output directory path (default: current directory)')
@click.option('--threads', '-t',
              default=4,
              type=int,
              callback=lambda ctx, param, value: _validate_threads(value) if value is not None else None,
              help='🧵 线程数 | Thread count (default: 4)')
@click.option('--mcyc-base',
              default=None,
              type=click.Path(),
              help='🗄️ MCycDB基础路径 | MCycDB base path')
@click.option('--skip-diamond',
              is_flag=True,
              help='⏭️ 跳过Diamond比对 | Skip Diamond alignment')
@click.option('--keep-temp',
              is_flag=True,
              help='💾 保留临时文件 | Keep temporary files')
def mcyc(input_list, output, threads, mcyc_base, skip_diamond, keep_temp):
    """
    甲烷循环基因丰度分析工具.

    基于MCycDB数据库分析甲烷循环相关基因的丰度，支持原始计数、TPM和CLR三种标准化方法。
    适用于宏基因组数据中甲烷循环功能基因的定量分析。

    输入文件格式 | Input File Format:
    样本名称和FastQ文件路径，用制表符分隔：

    sample1    /path/to/sample1_R1.fastq.gz
    sample2    /path/to/sample2_R1.fastq.gz;/path/to/sample2_R2.fastq.gz

    示例 | Examples:

    \b
    # 🎯 基本分析
    biopytools mcyc samples.txt

    \b
    # 🔧 指定输出目录和线程数
    biopytools mcyc samples.txt -o ./mcyc_results -t 8

    \b
    # ⚡ 跳过Diamond比对（如果已有中间结果）
    biopytools mcyc samples.txt --skip-diamond

    \b
    # 💾 保留临时文件用于调试
    biopytools mcyc samples.txt --keep-temp

    \b
    # 🗄️ 指定自定义MCycDB路径
    biopytools mcyc samples.txt --mcyc-base /path/to/MCycDB

    输出文件 | Output Files:
    - Matrix_1_RawCounts.txt: 原始计数矩阵
    - Matrix_2_TPM.txt: TPM标准化相对丰度矩阵
    - Matrix_3_CLR.txt: CLR变换矩阵（适用于GWAS分析）

    依赖工具 | Dependencies:
    - diamond: 序列比对工具
    - seqkit: 序列统计工具
    - perl: Perl解释器（用于运行MCycDB脚本）
    """

    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    mcyc_main = _lazy_import_mcyc_main()

    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['mcyc.py']

    # 必需参数 | Required parameter
    args.extend([input_list])

    # 可选参数 | Optional parameters
    if output is not None:
        args.extend(['--output', output])

    if threads != 4:
        args.extend(['--threads', str(threads)])

    if mcyc_base is not None:
        args.extend(['--mcyc-base', mcyc_base])

    if skip_diamond:
        args.append('--skip-diamond')

    if keep_temp:
        args.append('--keep-temp')

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 🚀 | Call original main function
        mcyc_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"❌ 错误 | Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv