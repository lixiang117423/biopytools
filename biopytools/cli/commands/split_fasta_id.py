"""
FASTA ID分割命令|FASTA ID Splitting Command
"""

import click
import sys
import os


def _lazy_import_split_fasta_id_main():
    """延迟加载split_fasta_id主函数|Lazy load split_fasta_id main function"""
    try:
        from ...split_fasta_id.main import main as split_fasta_id_main
        return split_fasta_id_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='FASTA ID分割工具|FASTA ID splitting tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTA文件路径|Input FASTA file path')
@click.option('--output', '-o',
              default='output.fasta',
              type=click.Path(),
              show_default=True,
              help='输出FASTA文件路径|Output FASTA file path')
@click.option('--position', '-p',
              default=0,
              type=int,
              show_default=True,
              help='提取位置(0表示第一个元素)|Extract position (0 means first element)')
@click.option('--delimiter', '-d',
              default='auto',
              type=str,
              show_default=True,
              help='分隔符类型|Delimiter type: "auto"(auto detect), "space", "tab", "both"(space and tab), or any character like "," or "|"')
@click.option('--keep-original',
              is_flag=True,
              help='保留原始文件作为备份|Keep original file as backup')
@click.option('--no-skip-empty',
              is_flag=True,
              help='不跳过空的序列名称行|Do not skip empty sequence name lines')
@click.option('--preserve-comments',
              is_flag=True,
              help='保留序列名称行中的注释|Preserve comments in sequence name lines')
def split_fasta_id(input, output, position, delimiter, keep_original,
                   no_skip_empty, preserve_comments):
    """
    FASTA ID分割工具|FASTA ID Splitting Tool

    从FASTA序列名称中提取指定位置的元素，简化序列ID|Extract element at specified position from FASTA sequence headers to simplify sequence IDs

    示例|Examples: biopytools split-fasta-id -i input.fasta -o output.fasta -p 0
    """

    # 延迟加载|Lazy load
    split_fasta_id_main = _lazy_import_split_fasta_id_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['split_fasta_id.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if output != 'output.fasta':
        args.extend(['-o', output])

    if position != 0:
        args.extend(['-p', str(position)])

    if delimiter != 'auto':
        args.extend(['-d', delimiter])

    # 处理选项|Processing options
    if keep_original:
        args.append('--keep-original')

    if no_skip_empty:
        args.append('--no-skip-empty')

    if preserve_comments:
        args.append('--preserve-comments')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        split_fasta_id_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
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
