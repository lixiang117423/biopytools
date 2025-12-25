"""
Rename Chromosomes Command
染色体重命名命令
"""

import click
import sys
import os


def _lazy_import_rename_main():
    """懒加载rename_chromosomes main函数 | Lazy load rename_chromosomes main function"""
    try:
        from ...rename_chromosomes.main import main as rename_main
        return rename_main
    except ImportError as e:
        click.echo(f"[ERROR] 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path


@click.command(short_help="染色体重命名工具",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='[FILE] 输入FASTA文件路径 | Input FASTA file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='[FILE] 输出FASTA文件路径 | Output FASTA file path')
@click.option('--number', '-n',
              required=True,
              type=int,
              help='[INT] 染色体数量 | Number of chromosomes')
def rename_chromosomes(input, output, number):
    """
    染色体重命名工具

    将FASTA文件的序列重命名，前N条序列命名为Chr01, Chr02...
    剩余序列命名为HiC_scaffold_01, HiC_scaffold_02...

    示例 | Examples:

    基本用法:
      biopytools rename-chromosomes -i genome.fa -o renamed.fa -n 20

    指定输出目录:
      biopytools rename-chromosomes -i input.fa -o output/renamed.fa -n 24

    输出说明 | Output Description:

    工具会生成以下文件:
      - renamed.fa: 重命名后的FASTA文件
      - rename_chromosomes_*.log: 运行日志

    命名规则 | Naming Rules:

      - 前 N 条序列命名为 Chr01, Chr02, ... ChrNN
      - 剩余序列命名为 HiC_scaffold_01, HiC_scaffold_02, ...

    注意事项 | Notes:

      - 染色体编号使用两位数格式 (01, 02, ..., 99)
      - Scaffold编号同样使用两位数格式
      - 输入文件必须是有效的FASTA格式
    """

    # 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    rename_main = _lazy_import_rename_main()

    # 构建参数列表传递给原始main函数 | Build argument list for original main function
    args = ['rename_chromosomes.py']

    # 必需参数 | Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])
    args.extend(['-n', str(number)])

    # 保存并恢复sys.argv | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 | Call original main function
        rename_main()
    except SystemExit as e:
        # 处理程序正常退出 | Handle normal program exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n[WARNING] 用户中断操作 | User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"[ERROR] 运行错误 | Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
