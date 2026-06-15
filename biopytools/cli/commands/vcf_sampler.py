"""
VCF抽样命令|VCF SNP Sampling Command
"""

import click
import sys
import os


def _lazy_import_vcf_sampler_main():
    """延迟加载vcf_sampler主函数|Lazy load vcf_sampler main function"""
    try:
        from ...vcf_sampler.main import main as vcf_sampler_main
        return vcf_sampler_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
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


def _validate_sample_rate(ctx, param, value):
    """验证抽样比例|Validate sample rate"""
    if value is not None and not (0.0 < value <= 1.0):
        raise click.BadParameter(f"抽样比例必须在0.0到1.0之间|Sample rate must be between 0.0 and 1.0, got: {value}")
    return value


@click.command(short_help="VCF抽样工具|VCF sampling tool",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入VCF文件路径|Input VCF file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出VCF文件路径|Output VCF file path')
@click.option('--sample-rate', '-r',
              type=float,
              default=0.25,
              show_default=True,
              callback=_validate_sample_rate,
              help='抽样比例|Sampling rate')
@click.option('--random-seed', '-s',
              type=int,
              default=1288,
              show_default=True,
              help='随机种子|Random seed')
@click.option('--log-file',
              type=click.Path(),
              default=None,
              help='日志文件路径|Log file path')
@click.option('--verbose', '-v',
              count=True,
              help='详细输出模式|Verbose mode')
def vcf_sampler(input, output, sample_rate, random_seed, log_file, verbose):
    """
    VCF抽样工具|VCF Sampling Tool

    从VCF文件中按比例随机抽取SNP位点|Randomly sample SNP sites from VCF files by proportion

    示例|Examples: biopytools vcf-sampler -i input.vcf.gz -o output.vcf.gz
    """

    # 延迟加载|Lazy loading
    vcf_sampler_main = _lazy_import_vcf_sampler_main()

    # 构建参数列表|Build argument list
    args = ['vcf_sampler.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数（只在非默认值时添加）|Optional parameters (add only when non-default)
    if sample_rate != 0.25:
        args.extend(['-r', str(sample_rate)])

    if random_seed != 1288:
        args.extend(['-s', str(random_seed)])
    else:
        # 总是添加随机种子参数|Always add random seed parameter
        args.extend(['-s', str(random_seed)])

    if log_file:
        args.extend(['--log-file', log_file])

    # 详细输出模式|Verbose mode
    if verbose == 1:
        args.append('-v')
    elif verbose >= 2:
        args.append('-vv')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        vcf_sampler_main()
    except SystemExit as e:
        # 处理程序正常退出|Handle normal program exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(130)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
