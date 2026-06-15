"""
EGAPx批量运行配置生成命令|EGAPx Batch Config Generator Command
"""

import click
import sys
import os
from ...common.paths import expand_path


def _lazy_import_egapx_batch_main():
    """延迟加载egapx_batch主函数|Lazy load egapx-batch main function"""
    try:
        from ...egapx_batch.main import main as egapx_batch_main
        return egapx_batch_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_genome(file_path):
    """验证基因组文件(仅在非帮助模式)|Validate genome file (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not file_path:
        raise click.BadParameter("基因组文件路径不能为空|Genome file path cannot be empty")

    # 展开路径|Expand path
    file_path = expand_path(file_path)

    if not os.path.exists(file_path):
        raise click.BadParameter(f"基因组文件不存在|Genome file not found: {file_path}")
    if not file_path.endswith(('.fa', '.fa.gz', '.fasta', '.fasta.gz')):
        raise click.BadParameter(f"无效的基因组文件格式|Invalid genome file format: {file_path}")
    return file_path


def _validate_egapx_path(egapx_path):
    """验证EGAPx路径(仅在非帮助模式)|Validate EGAPx path (only in non-help mode)"""
    if _is_help_request():
        return egapx_path
    if not egapx_path:
        raise click.BadParameter("EGAPx路径不能为空|EGAPx path cannot be empty")

    # 展开路径|Expand path
    egapx_path = expand_path(egapx_path)

    if not os.path.exists(egapx_path):
        raise click.BadParameter(f"EGAPx路径不存在|EGAPx path not found: {egapx_path}")
    # EGAPx可以是一个目录或一个可执行文件|EGAPx can be a directory or an executable file
    if not (os.path.isdir(egapx_path) or os.path.isfile(egapx_path)):
        raise click.BadParameter(f"EGAPx路径无效|EGAPx path is invalid: {egapx_path}")
    return egapx_path


def _validate_dir_path(dir_path):
    """验证目录路径(仅在非帮助模式)|Validate directory path (only in non-help mode)"""
    if _is_help_request():
        return dir_path
    if not dir_path:
        return dir_path  # 空值是允许的（可选参数）|Empty value is allowed (optional parameter)

    # 展开路径|Expand path
    dir_path = expand_path(dir_path)

    if not os.path.exists(dir_path):
        raise click.BadParameter(f"目录不存在|Directory not found: {dir_path}")
    if not os.path.isdir(dir_path):
        raise click.BadParameter(f"路径不是目录|Path is not a directory: {dir_path}")
    return dir_path


def _validate_reads_file(file_path):
    """验证测序数据文件(仅在非帮助模式)|Validate reads file (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not file_path:
        return file_path  # 空值是允许的（可选参数）|Empty value is allowed (optional parameter)

    # 展开路径|Expand path
    file_path = expand_path(file_path)

    if not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File not found: {file_path}")
    # 支持fq、fq.gz、fastq、fastq.gz格式|Support fq, fq.gz, fastq, fastq.gz formats
    if not file_path.endswith(('.fq', '.fq.gz', '.fastq', '.fastq.gz')):
        raise click.BadParameter(f"无效的测序数据文件格式，应为.fq或.fastq格式|Invalid reads file format, expected .fq or .fastq: {file_path}")
    return file_path


def _validate_reads_input(input_path):
    """验证测序数据输入(文件或目录，仅在非帮助模式)|Validate reads input (file or directory, only in non-help mode)"""
    if _is_help_request():
        return input_path
    if not input_path:
        return input_path  # 空值是允许的（可选参数）|Empty value is allowed (optional parameter)

    # 展开路径|Expand path
    input_path = expand_path(input_path)

    if not os.path.exists(input_path):
        raise click.BadParameter(f"路径不存在|Path not found: {input_path}")

    # 可以是文件或目录|Can be file or directory
    if not (os.path.isfile(input_path) or os.path.isdir(input_path)):
        raise click.BadParameter(f"路径既不是文件也不是目录|Path is neither file nor directory: {input_path}")

    return input_path


def _validate_sif_image(file_path):
    """验证Singularity镜像文件(仅在非帮助模式)|Validate Singularity image file (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not file_path:
        return file_path

    file_path = expand_path(file_path)

    if not os.path.exists(file_path):
        raise click.BadParameter(f"SIF镜像文件不存在|SIF image not found: {file_path}")
    if not os.path.isfile(file_path):
        raise click.BadParameter(f"SIF路径不是文件|SIF path is not a file: {file_path}")

    return file_path


@click.command(
    short_help='EGAPx批量运行配置生成工具|EGAPx Batch Config Generator',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_genome(value) if value else None,
              help='基因组FASTA文件路径|Genome FASTA file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录路径|Output directory path')
@click.option('--egapx', '-e',
              default='~/software/EGAPX_v.0.4.1-alpha/egapx',
              callback=lambda ctx, param, value: _validate_egapx_path(value) if value else None,
              show_default=True,
              help='EGAPx安装路径|EGAPx installation path')
@click.option('--local-cache',
              default='~/software/EGAPX_v.0.4.1-alpha/local_cache',
              callback=lambda ctx, param, value: _validate_dir_path(value) if value else value,
              show_default=True,
              help='EGAPx本地缓存路径|EGAPx local cache path')
@click.option('--sif',
              default='~/software/EGAPX_v.0.4.1-alpha/egapx/egapx_0.4.1-alpha.sif',
              callback=lambda ctx, param, value: _validate_sif_image(value) if value else value,
              show_default=True,
              help='Singularity镜像路径|Singularity image path')
@click.option('--no-split',
              is_flag=True,
              default=False,
              show_default=True,
              help='不按染色体拆分基因组|Do not split genome by chromosome')
@click.option('--chr-prefix', '-p',
              help='染色体前缀过滤|Chromosome prefix filter')
@click.option('--locus-prefix',
              default='',
              show_default=True,
              help='locus标签前缀|Locus tag prefix')
@click.option('--report-name',
              default='EGAPx',
              show_default=True,
              help='报告名称|Report name')
@click.option('--short-reads',
              default='',
              callback=lambda ctx, param, value: _validate_reads_input(value) if value else value,
              show_default=True,
              help='短读长测序数据(目录或文件)|Short reads (directory or file)')
@click.option('--long-reads',
              default='',
              callback=lambda ctx, param, value: _validate_reads_input(value) if value else value,
              show_default=True,
              help='长读长测序数据(目录或文件)|Long reads (directory or file)')
@click.option('--taxid',
              default='71234',
              show_default=True,
              help='物种分类ID|Species taxonomy ID')
def egapx_batch(genome, output, egapx, local_cache, sif, no_split, chr_prefix, locus_prefix,
                report_name, short_reads, long_reads, taxid):
    """
    EGAPx批量运行配置生成工具|EGAPx Batch Config Generator

    按染色体拆分基因组并批量生成EGAPx运行配置和脚本|Split genome by chromosome and generate EGAPx configs and scripts

    示例|Examples: biopytools egapx-batch -g genome.fa -o output_dir
    """

    # 延迟加载|Lazy load: import only when actually called
    egapx_batch_main = _lazy_import_egapx_batch_main()

    # 构建主函数的参数列表|Build argument list for main function
    args = ['egapx_batch.py']

    # 必需参数|Required parameters
    args.extend(['-g', genome])
    args.extend(['-o', output])

    # 可选参数|Optional parameters (add only when non-default)
    if egapx != '~/software/EGAPX_v.0.4.1-alpha/egapx':
        args.extend(['-e', egapx])

    if local_cache != '~/software/EGAPX_v.0.4.1-alpha/local_cache':
        args.extend(['--local-cache', local_cache])

    if sif:
        args.extend(['--sif', sif])

    if no_split:
        args.append('--no-split')

    if chr_prefix:
        args.extend(['-p', chr_prefix])

    if locus_prefix:
        args.extend(['--locus-prefix', locus_prefix])

    if report_name != 'EGAPx':
        args.extend(['--report-name', report_name])

    if short_reads:
        args.extend(['--short-reads', short_reads])

    if long_reads:
        args.extend(['--long-reads', long_reads])

    if taxid != '71234':
        args.extend(['--taxid', taxid])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        egapx_batch_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
