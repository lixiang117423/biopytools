"""
InterProScan蛋白质功能注释|InterProScan Protein Function Annotation Command
"""

import click
import sys
import os


def _lazy_import_interproscan_main():
    """延迟加载interproscan主函数|Lazy load interproscan main function"""
    try:
        from ...interproscan.main import main as interproscan_main
        return interproscan_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='InterProScan蛋白质功能注释工具|InterProScan protein function annotation tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='蛋白质FASTA文件|Protein FASTA file path')
@click.option('--output-prefix', '-o',
              required=True,
              help='输出文件前缀|Output file prefix (without extension)')
@click.option('--interproscan-path', '-a',
              default='~/software/InterProScan/v.5.75-106.0/interproscan-5.75-106.0/interproscan.sh',
              show_default=True,
              help='InterProScan软件路径|InterProScan software installation path')
@click.option('--format', '-f',
              default='TSV',
              type=click.Choice(['TSV', 'GFF3', 'XML', 'HTML', 'JSON', 'TXT']),
              show_default=True,
              help='输出格式|Output format')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--disable-precalc',
              is_flag=True,
              default=True,
              help='禁用预计算查找服务(默认)|Disable precalc lookup service (default)')
@click.option('--enable-precalc',
              is_flag=True,
              help='启用预计算查找服务(需要网络)|Enable precalc lookup service (requires network)')
@click.option('--goterms',
              is_flag=True,
              help='获取GO术语(默认不获取)|Get GO terms (disabled by default)')
@click.option('--no-goterms',
              is_flag=True,
              default=True,
              help='不获取GO术语(默认)|Do not get GO terms (default)')
@click.option('--pathways',
              is_flag=True,
              help='获取Pathway信息(默认不获取)|Get pathway information (disabled by default)')
@click.option('--no-pathways',
              is_flag=True,
              default=True,
              help='不获取Pathway信息(默认)|Do not get pathway information (default)')
@click.option('--applications', '-appl',
              help='指定运行的应用，逗号分隔|Specify applications to run, comma-separated')
@click.option('--temp-dir',
              type=click.Path(),
              help='临时目录路径|Temporary directory path')
@click.option('--python-path',
              type=click.Path(),
              help='Python解释器路径（兼容Python 3.8-3.11）|Python interpreter path (compatible with Python 3.8-3.11)')
def interproscan(input, output_prefix, interproscan_path, format, threads,
                disable_precalc, enable_precalc, goterms, no_goterms,
                pathways, no_pathways, applications, temp_dir, python_path):
    """
    InterProScan蛋白质功能注释工具|InterProScan Protein Function Annotation Tool

    对蛋白质序列进行功能结构域注释和GO术语预测|Functional domain annotation and GO term prediction for protein sequences

    示例|Example: biopytools interproscan -i proteins.fa -o results -t 32
    """

    # 延迟加载|Lazy loading
    interproscan_main = _lazy_import_interproscan_main()

    # 构建参数列表|Build argument list
    args = ['interproscan.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output_prefix])

    # 可选参数|Optional parameters
    if interproscan_path != '~/software/InterProScan/v.5.75-106.0/interproscan-5.75-106.0/interproscan.sh':
        args.extend(['--interproscan-path', interproscan_path])

    if format != 'TSV':
        args.extend(['--format', format])

    if threads != 24:
        args.extend(['--threads', str(threads)])

    # 处理选项|Processing options
    if enable_precalc:
        args.append('--enable-precalc')

    if goterms:
        args.append('--goterms')

    if pathways:
        args.append('--pathways')

    if applications:
        args.extend(['--applications', applications])

    if temp_dir:
        args.extend(['--temp-dir', temp_dir])

    if python_path:
        args.extend(['--python-path', python_path])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        interproscan_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
