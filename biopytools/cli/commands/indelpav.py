"""
INDEL PAV分析CLI包装器|INDEL PAV Analysis CLI Wrapper
"""

import click
import sys
import os


def _lazy_import_pav_main():
    """延迟加载PAV主函数|Lazy load PAV main function"""
    try:
        from ...indelpav.main import main as pav_main
        return pav_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(short_help='INDEL PAV分析|INDEL PAV Analysis',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--vcf', '-v',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入VCF文件路径|Input VCF file path')
@click.option('--output', '-o',
              default='./indel_pav.txt',
              type=click.Path(),
              show_default=True,
              help='输出文件路径|Output file path')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--min-length',
              type=int,
              default=1,
              show_default=True,
              help='最小INDEL长度(bp)|Minimum INDEL length (bp)')
@click.option('--max-length',
              type=int,
              help='最大INDEL长度(bp)|Maximum INDEL length (bp)')
@click.option('--min-quality', '-q',
              type=float,
              default=20.0,
              show_default=True,
              help='最小质量分数|Minimum quality score')
@click.option('--min-depth', '-d',
              type=int,
              default=5,
              show_default=True,
              help='最小深度|Minimum depth')
@click.option('--max-missing',
              type=float,
              default=0.8,
              show_default=True,
              help='最大缺失率(0-1)|Maximum missing rate (0-1)')
@click.option('--include-complex',
              is_flag=True,
              help='包含复杂变异|Include complex variants')
@click.option('--compress',
              is_flag=True,
              help='压缩输出文件|Compress output file')
@click.option('--bcftools-path',
              default='bcftools',
              show_default=True,
              help='BCFtools软件路径|BCFtools software path')
def indelpav(vcf, output, threads, min_length, max_length, min_quality, min_depth,
              max_missing, include_complex, compress, bcftools_path):
    """
    INDEL PAV分析流程|INDEL PAV Analysis Pipeline

    基于VCF文件的INDEL存在/缺失变异分析
    INDEL presence/absence variation analysis based on VCF file

    示例|Examples: biopytools indelpav -v variants.vcf -o pav_results.txt
    """

    # 延迟加载|Lazy load
    pav_main = _lazy_import_pav_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['indelpav.py']

    # 必需参数|Required parameters
    args.extend(['-v', vcf])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if output != './indel_pav.txt':
        args.extend(['-o', output])

    if threads != 88:
        args.extend(['-t', str(threads)])

    if min_length != 1:
        args.extend(['--min-length', str(min_length)])

    if max_length is not None:
        args.extend(['--max-length', str(max_length)])

    if min_quality != 20.0:
        args.extend(['-q', str(min_quality)])

    if min_depth != 5:
        args.extend(['-d', str(min_depth)])

    if max_missing != 0.8:
        args.extend(['--max-missing', str(max_missing)])

    if bcftools_path != 'bcftools':
        args.extend(['--bcftools-path', bcftools_path])

    # 布尔选项|Boolean options
    if include_complex:
        args.append('--include-complex')

    if compress:
        args.append('--compress')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        pav_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"意外错误|Unexpected error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
