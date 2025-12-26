"""
Assembly to AGP Converter Command
Assembly文件转换为AGP格式并生成染色体列表文件
"""

import click
import sys
import os


def _lazy_import_assembly2agp_main():
    """懒加载assembly2agp main函数 | Lazy load assembly2agp main function"""
    try:
        from ...assembly2agp.main import main as assembly2agp_main
        return assembly2agp_main
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


@click.command(short_help="Assembly文件转AGP格式工具",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--assembly', '-a',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📂 Assembly文件路径 | Assembly file path')
@click.option('--prefix', '-p',
              required=True,
              help='🏷️ 输出文件前缀 | Output prefix (for both AGP and chr.list files)')
@click.option('--output-dir', '-o',
              default='.',
              type=click.Path(),
              help='📁 输出目录路径 | Output directory path (default: current directory)')
@click.option('--gap', '-g',
              type=int,
              default=100,
              help='📏 Scaffold间gap大小 (bp) | Gap size between scaffolds in bp (default: 100)')
@click.option('--num-chromosomes', '-n',
              required=True,
              type=int,
              help='🧬 染色体数量 | Number of chromosomes to output')
@click.option('--force', '-f',
              is_flag=True,
              help='🔄 强制覆盖已存在的输出文件 | Force overwrite existing output files')
@click.option('--verbose', '-v',
              count=True,
              help='📢 详细输出模式 (-v: INFO, -vv: DEBUG) | Verbose mode (-v: INFO, -vv: DEBUG)')
@click.option('--quiet',
              is_flag=True,
              help='🔇 静默模式 (只输出ERROR) | Quiet mode (only ERROR)')
@click.option('--log-file',
              type=click.Path(),
              help='📝 日志文件路径 | Log file path')
def assembly2agp(assembly, prefix, output_dir, gap, num_chromosomes,
                 force, verbose, quiet, log_file):
    """
    Assembly文件转AGP格式工具

    将assembly文件转换为AGP格式，并生成染色体列表文件。

    示例 | Examples:

    \b
    # 基本使用
    biopytools assembly2agp \\
        -a corrected_asm.FINAL.assembly \\
        -p output_prefix \\
        -n 12

    \b
    # 指定输出目录和gap大小
    biopytools assembly2agp \\
        -a corrected_asm.FINAL.assembly \\
        -p output_prefix \\
        -o ./output \\
        -g 200 \\
        -n 10

    \b
    # 强制覆盖已存在的文件
    biopytools assembly2agp \\
        -a corrected_asm.FINAL.assembly \\
        -p output_prefix \\
        -n 12 \\
        -f

    \b
    # 启用详细输出
    biopytools assembly2agp \\
        -a corrected_asm.FINAL.assembly \\
        -p output_prefix \\
        -n 12 \\
        -v
    """

    # 验证染色体数量
    if num_chromosomes <= 0:
        click.echo("❌ 错误 | Error: 染色体数量必须大于0 | Number of chromosomes must be greater than 0", err=True)
        sys.exit(1)

    # 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    assembly2agp_main = _lazy_import_assembly2agp_main()

    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['assembly2agp.py']

    # 必需参数 📋 | Required parameters
    args.extend(['-a', assembly])
    args.extend(['-p', prefix])

    # 可选参数（只在非默认值时添加，减少命令行长度）⚙️ | Optional parameters (add only when non-default)
    if output_dir != '.':
        args.extend(['-o', output_dir])

    if gap != 100:
        args.extend(['-g', str(gap)])

    args.extend(['-n', str(num_chromosomes)])

    # 处理选项（布尔标志）🚩 | Processing options (boolean flags)
    if force:
        args.append('--force')

    # 日志参数 📝 | Logging parameters
    if verbose:
        args.extend(['-v'] * verbose)

    if quiet:
        args.append('--quiet')

    if log_file:
        args.extend(['--log-file', log_file])

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 🚀 | Call original main function
        assembly2agp_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"❌ 错误 | Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
