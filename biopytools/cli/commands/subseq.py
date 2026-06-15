"""
序列子集提取命令|Sequence Subsequence Extraction Command
"""

import click
import sys
import os


def _lazy_import_subseq_main():
    """延迟加载subseq主函数|Lazy load subseq main function"""
    try:
        from ...subseq.main import main as subseq_main
        return subseq_main
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
    short_help='序列子集提取工具|Sequence subsequence extraction tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTA文件路径|Input FASTA file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出FASTA文件路径|Output FASTA file path')
@click.option('--id-list', '-l',
              type=click.Path(exists=True),
              help='ID列表文件路径|ID list file path')
@click.option('--pattern', '-p',
              help='模式匹配字符串|Pattern matching string')
@click.option('--pattern-type',
              type=click.Choice(['contains', 'startswith', 'endswith', 'regex']),
              default='contains',
              show_default=True,
              help='模式类型|Pattern type')
@click.option('--ignore-case',
              is_flag=True,
              help='忽略大小写|Case insensitive')
@click.option('--min-length',
              type=int,
              default=0,
              show_default=True,
              help='最小序列长度|Minimum sequence length')
@click.option('--max-length',
              type=int,
              help='最大序列长度|Maximum sequence length')
@click.option('--no-order',
              is_flag=True,
              help='不保持ID列表顺序|Do not keep ID list order')
@click.option('--log-dir',
              default='.',
              type=click.Path(),
              show_default=True,
              help='日志输出目录|Log output directory')
def subseq(input, output, id_list, pattern, pattern_type, ignore_case,
           min_length, max_length, no_order, log_dir):
    """
    序列子集提取工具|Sequence Subsequence Extraction Tool

    根据ID列表、模式匹配或长度范围从FASTA文件中提取序列子集|Extract sequence subsets from FASTA files by ID list, pattern matching, or length range

    示例|Examples: biopytools subseq -i input.fasta -l id_list.txt -o output.fasta
    """

    # 延迟加载|Lazy load
    subseq_main = _lazy_import_subseq_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['subseq.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 提取方式检查|Extraction method check
    extraction_methods = 0
    if id_list:
        extraction_methods += 1
    if pattern:
        extraction_methods += 1
    if not id_list and not pattern:
        # 使用长度筛选|Use length filtering
        extraction_methods += 1

    if extraction_methods > 1:
        click.echo("错误|Error: Only one extraction method can be specified", err=True)
        sys.exit(1)

    # 根据提取方式添加参数|Add parameters based on extraction method
    if id_list:
        args.extend(['-l', id_list])
    elif pattern:
        args.extend(['-p', pattern])
        if pattern_type != 'contains':
            args.extend(['--pattern-type', pattern_type])
        if ignore_case:
            args.append('--ignore-case')
    else:
        # 长度筛选|Length filtering
        args.append('--length-only')
        if min_length != 0:
            args.extend(['--min-length', str(min_length)])
        if max_length:
            args.extend(['--max-length', str(max_length)])

    # 其他选项|Other options
    if no_order:
        args.append('--no-order')

    if log_dir != '.':
        args.extend(['--log-dir', log_dir])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        subseq_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
