"""
 BAM深度分析命令|BAM Depth Analysis Command
"""

import click
import sys
import os


def _lazy_import_depth_main():
    """延迟加载depth主函数|Lazy load depth main function"""
    try:
        from ...coverage.main import main as depth_main
        return depth_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_multiple_files(files_tuple):
    """验证多个文件存在(仅在非帮助模式)|Validate multiple files existence (only in non-help mode)"""
    if _is_help_request():
        return files_tuple

    for file_path in files_tuple:
        if not os.path.exists(file_path):
            raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")

    return files_tuple


@click.command(short_help=" BAM/SAM深度分析",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              multiple=True,
              required=True,
              callback=lambda ctx, param, value: _validate_multiple_files(value) if value else None,
              help='BAM/SAM文件|Input BAM/SAM file paths or directories')
@click.option('--output', '-o',
              default='./depth_results.txt',
              type=click.Path(),
              help='输出文件|Output file path or directory (default: ./depth_results.txt)')
@click.option('--chromosome', '-c',
              default='all',
              help='目标染色体|Target chromosome name (default: all)')
@click.option('--region', '-r',
              default='all',
              help='染色体区域:start:end|Chromosome region, format: start:end (default: all)')
@click.option('--threads', '-t',
              type=int,
              default=12,
              help='线程数|Number of threads (default: 88)')
@click.option('--quality', '-q',
              type=int,
              default=0,
              help='最小碱基质量|Minimum base quality threshold (default: 0)')
@click.option('--mapping-quality', '-Q',
              type=int,
              default=0,
              help='最小比对质量|Minimum mapping quality threshold (default: 0)')
@click.option('--samtools-path',
              default='samtools',
              help='samtools路径|samtools software path (default: samtools)')
@click.option('--output-format',
              type=click.Choice(['txt']),
              default='txt',
              help='输出格式|Output format (default: txt)')
@click.option('--compress',
              is_flag=True,
              help='压缩输出文件|Compress output file')
@click.option('--enable-windows',
              is_flag=True,
              help='启用滑动窗口分析|Enable sliding window analysis')
@click.option('--window-size',
              type=int,
              default=1000,
              help='窗口大小(bp)|Window size (bp) (default: 1000)')
@click.option('--window-step',
              type=int,
              default=0,
              help='窗口步长(bp)|Window step size (bp) (default: 0)')
def coverage(input, output, chromosome, region, threads, quality, mapping_quality,
             samtools_path, output_format, compress, enable_windows, window_size, window_step):
    """
     BAM/SAM深度分析工具

    计算BAM/SAM文件的测序深度
    支持多文件并行处理和滑动窗口分析

    使用示例|Examples:

    \b
    # 基本用法 - 分析BAM文件深度
    biopytools coverage -i sample.bam -o depth_results.txt

    \b
    # 指定染色体区域
    biopytools coverage -i sample.bam -o results.txt \\
        -c chr12 -r 136491092:138554123

    \b
    # 多文件并行处理
    biopytools coverage -i sample1.bam sample2.bam \\
        -o results.txt --threads 16

    \b
    # 滑动窗口分析
    biopytools coverage -i sample.bam -o results.txt \\
        --enable-windows --window-size 1000 --window-step 500

    \b
    # 质量过滤
    biopytools coverage -i data.bam -o results.txt \\
        -q 20 -Q 30 --compress
    """

    # 延迟加载: 仅在实际调用时导入
    depth_main = _lazy_import_depth_main()

    # 构建主函数的参数列表
    args = ['coverage.py']

    # 必需参数
    args.extend(['-i'] + list(input))

    # 可选参数
    if output != './depth_results.txt':
        args.extend(['-o', output])

    if chromosome != 'all':
        args.extend(['-c', chromosome])

    if region != 'all':
        args.extend(['-r', region])

    if threads != 88:
        args.extend(['-t', str(threads)])

    if quality != 0:
        args.extend(['-q', str(quality)])

    if mapping_quality != 0:
        args.extend(['-Q', str(mapping_quality)])

    if samtools_path != 'samtools':
        args.extend(['--samtools-path', samtools_path])

    if output_format != 'txt':
        args.extend(['--output-format', output_format])

    # 布尔选项
    if compress:
        args.append('--compress')

    if enable_windows:
        args.append('--enable-windows')

    if window_size != 1000:
        args.extend(['--window-size', str(window_size)])

    if window_step != 0:
        args.extend(['--window-step', str(window_step)])

    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数
        depth_main()
    except SystemExit as e:
        # 处理正常程序退出
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv