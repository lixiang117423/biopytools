"""
VCF合并命令|VCF Merge by Chromosome Command
"""

import click
import sys
import os


def _lazy_import_vcf_merger_main():
    """延迟加载vcf_merger主函数|Lazy load vcf_merger main function"""
    try:
        from ...vcf_merger.main import main as vcf_merger_main
        return vcf_merger_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_dir_exists(dir_path):
    """验证目录存在(仅在非帮助模式)|Validate directory existence (only in non-help mode)"""
    if not _is_help_request() and dir_path and not os.path.exists(dir_path):
        raise click.BadParameter(f"目录不存在|Directory does not exist: {dir_path}")
    return dir_path


@click.command(
    short_help='VCF按染色体合并工具|VCF merge by chromosome tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input-dir', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_dir_exists(value) if value else None,
              help='输入VCF文件目录|Input directory containing VCF files')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--pattern',
              default='*.joint.vcf.gz',
              show_default=True,
              help='VCF文件名模式|VCF file name pattern')
@click.option('--no-index',
              is_flag=True,
              help='不生成索引文件|Do not generate index files')
@click.option('--verbose', '-v',
              count=True,
              help='详细输出模式|Verbose mode')
@click.option('--quiet',
              is_flag=True,
              help='静默模式|Quiet mode')
def vcf_merger(input_dir, output_dir, threads, pattern, no_index, verbose, quiet):
    """
    VCF按染色体合并工具|VCF Merge by Chromosome Tool

    使用bcftools按染色体合并VCF文件|Merge VCF files by chromosome using bcftools

    示例|Examples: biopytools vcf-merger -i /path/to/vcf_files -o /path/to/output
    """

    # 延迟加载|Lazy load
    vcf_merger_main = _lazy_import_vcf_merger_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['vcf_merger.py']

    # 必需参数|Required parameters
    args.extend(['-i', input_dir])
    args.extend(['-o', output_dir])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if threads != 4:
        args.extend(['-t', str(threads)])

    if pattern != '*.joint.vcf.gz':
        args.extend(['--pattern', pattern])

    if no_index:
        args.append('--no-index')

    # 日志参数|Logging parameters
    if verbose > 0:
        args.extend(['-v'] * verbose)

    if quiet:
        args.append('--quiet')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        vcf_merger_main()
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
