"""
Assembly转AGP格式命令|Assembly to AGP Converter Command
"""

import click
import sys
import os


def _lazy_import_assembly2agp_main():
    """延迟加载assembly2agp主函数|Lazy load assembly2agp main function"""
    try:
        from ...assembly2agp.main import main as assembly2agp_main
        return assembly2agp_main
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
    short_help='Assembly转AGP格式|Assembly to AGP format converter',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--assembly', '-a',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='组装文件|Assembly file path')
@click.option('--prefix', '-p',
              required=True,
              help='输出前缀|Output prefix')
@click.option('--output-dir', '-o',
              default='.',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory path')
@click.option('--gap', '-g',
              type=int,
              default=100,
              show_default=True,
              help='Scaffold间隙(bp)|Scaffold gap size in bp')
@click.option('--num-chromosomes', '-n',
              required=True,
              type=int,
              help='染色体数量|Number of chromosomes')
@click.option('--force', '-f',
              is_flag=True,
              help='强制覆盖|Force overwrite existing files')
@click.option('--verbose', '-v',
              count=True,
              help='详细模式(-v: INFO, -vv: DEBUG)|Verbose mode (-v: INFO, -vv: DEBUG)')
@click.option('--quiet',
              is_flag=True,
              help='静默模式(仅ERROR)|Quiet mode (ERROR only)')
@click.option('--log-file',
              type=click.Path(),
              help='日志文件|Log file path')
def assembly2agp(assembly, prefix, output_dir, gap, num_chromosomes,
                 force, verbose, quiet, log_file):
    """
    Assembly转AGP格式工具|Assembly to AGP Format Converter

    将Assembly格式转换为AGP格式并生成染色体列表|Convert assembly format to AGP format and generate chromosome list

    示例|Examples: biopytools assembly2agp -a assembly.fasta -p output_prefix -n 12
    """

    # 参数验证|Parameter validation
    if num_chromosomes <= 0:
        click.echo("错误|Error: 染色体数量必须大于0|Number of chromosomes must be greater than 0", err=True)
        sys.exit(1)

    # 延迟加载|Lazy loading
    assembly2agp_main = _lazy_import_assembly2agp_main()

    # 构建参数列表|Build argument list
    args = ['assembly2agp.py']

    # 必需参数|Required parameters
    args.extend(['-a', assembly])
    args.extend(['-p', prefix])

    # 可选参数|Optional parameters
    if output_dir != '.':
        args.extend(['-o', output_dir])

    if gap != 100:
        args.extend(['-g', str(gap)])

    args.extend(['-n', str(num_chromosomes)])

    # 布尔选项|Boolean options
    if force:
        args.append('--force')

    # 日志参数|Logging parameters
    if verbose:
        args.extend(['-v'] * verbose)

    if quiet:
        args.append('--quiet')

    if log_file:
        args.extend(['--log-file', log_file])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        assembly2agp_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
