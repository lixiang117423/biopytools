"""
🧬 序列子集提取命令 | Sequence Subsequence Extraction Command
🚀 高级优化版本：解决--help响应速度问题 ⚡️
"""

import click
import sys
import os


def _lazy_import_subseq_main():
    """😴 懒加载subseq main函数 | Lazy load subseq main function"""
    try:
        from ...subseq.main import main as subseq_main
        return subseq_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """❓ 检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """🔍 验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"📂❌ 文件不存在 | File does not exist: {file_path}")
    return file_path


@click.command(short_help="🧬 序列子集提取",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📄 输入FASTA文件路径 | Input FASTA file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='📤 输出FASTA文件路径 | Output FASTA file path')
@click.option('--id-list', '-l',
              type=click.Path(exists=True),
              help='📋 ID列表文件路径 | ID list file path')
@click.option('--pattern', '-p',
              help='🔍 模式匹配字符串 | Pattern matching string')
@click.option('--pattern-type',
              type=click.Choice(['contains', 'startswith', 'endswith', 'regex']),
              default='contains',
              help='📝 模式类型 | Pattern type (default: contains)')
@click.option('--ignore-case',
              is_flag=True,
              help='🔤 忽略大小写 | Case insensitive')
@click.option('--min-length',
              type=int,
              default=0,
              help='📏 最小序列长度 | Minimum sequence length (default: 0)')
@click.option('--max-length',
              type=int,
              help='📏 最大序列长度 | Maximum sequence length (default: unlimited)')
@click.option('--no-order',
              is_flag=True,
              help='🔄 不保持ID列表顺序 | Do not keep ID list order')
@click.option('--log-dir',
              default='.',
              type=click.Path(),
              help='📁 日志输出目录 | Log output directory (default: current directory)')
def subseq(input, output, id_list, pattern, pattern_type, ignore_case,
           min_length, max_length, no_order, log_dir):
    """
    🧬 序列子集提取工具.

    💡 支持三种提取方式 | Supports three extraction methods:

    1. 📋 根据ID列表提取 | Extract by ID list
    2. 🔍 根据模式匹配提取 | Extract by pattern matching
    3. 📏 根据长度范围提取 | Extract by length range

    📚 示例 | Examples:

    \b
    # 📋 根据ID列表提取序列
    biopytools subseq -i input.fasta -l id_list.txt -o output.fasta

    \b
    # 🔍 根据模式匹配提取序列（以特定前缀开头）
    biopytools subseq -i input.fasta -p "chr1_" -o output.fasta --pattern-type startswith

    \b
    # 🔍 根据模式匹配提取序列（忽略大小写）
    biopytools subseq -i input.fasta -p "GENE" -o output.fasta --ignore-case

    \b
    # 📏 根据长度范围提取序列
    biopytools subseq -i input.fasta -o output.fasta --length-only --min-length 1000 --max-length 5000

    \b
    # 🔄 不保持ID列表顺序（使用FASTA文件中的原始顺序）
    biopytools subseq -i input.fasta -l id_list.txt -o output.fasta --no-order
    """

    # 😴 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    subseq_main = _lazy_import_subseq_main()

    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['subseq.py']

    # 必需参数 📋 | Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 提取方式检查 | Extraction method check
    extraction_methods = 0
    if id_list:
        extraction_methods += 1
    if pattern:
        extraction_methods += 1
    if not id_list and not pattern:
        # 使用长度筛选 | Use length filtering
        extraction_methods += 1

    if extraction_methods > 1:
        click.echo("❌ 只能指定一种提取方式 | Only one extraction method can be specified", err=True)
        sys.exit(1)

    # 根据提取方式添加参数 | Add parameters based on extraction method
    if id_list:
        args.extend(['-l', id_list])
    elif pattern:
        args.extend(['-p', pattern])
        args.extend(['--pattern-type', pattern_type])
        if ignore_case:
            args.append('--ignore-case')
    else:
        # 长度筛选 | Length filtering
        args.append('--length-only')
        args.extend(['--min-length', str(min_length)])
        if max_length:
            args.extend(['--max-length', str(max_length)])

    # 其他选项 | Other options
    if no_order:
        args.append('--no-order')

    if log_dir != '.':
        args.extend(['--log-dir', log_dir])

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 🚀 | Call original main function
        subseq_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"❌ 运行错误 | Runtime Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv