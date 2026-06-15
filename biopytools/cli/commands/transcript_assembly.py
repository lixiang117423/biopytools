"""
转录本从头组装命令|Transcript De Novo Assembly Command
"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...transcript_assembly.main import main as ta_main
        return ta_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path exists (only in non-help mode)"""
    if _is_help_request():
        return path
    if not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='转录本从头组装|Transcript de novo assembly',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='基因组FASTA文件路径|Genome FASTA file path')
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='输入FASTQ文件目录|Input FASTQ file directory')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--pattern', '-p',
              type=str,
              default='*_1.clean.fq.gz',
              show_default=True,
              help='FASTQ文件命名模式|FASTQ file naming pattern (* is sample name placeholder)')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--sample-timeout',
              type=int,
              default=43200,
              show_default=True,
              help='单个样本处理超时时间（秒）|Sample processing timeout in seconds')
@click.option('--step', '-s',
              type=click.Choice(['1', '2', '3', '4', '5', '6']),
              help='运行指定步骤|Run only specified step')
@click.option('--verbose', '-v',
              count=True,
              help='增加输出详细程度|Increase output verbosity')
@click.option('--quiet',
              is_flag=True,
              help='静默模式|Quiet mode, only output errors')
@click.option('--force',
              is_flag=True,
              help='强制重新处理已完成的步骤|Force re-process completed steps')
def transcript_assembly(genome, input, output, pattern, threads, sample_timeout,
                         step, verbose, quiet, force):
    """
    转录本从头组装流程：HISAT2 + StringTie|Transcript de novo assembly pipeline: HISAT2 + StringTie

    基于参考基因组的转录本从头组装，适用于无GFF注释的场景|De novo transcript assembly based on reference genome, suitable for scenarios without GFF annotation

    示例|Examples: biopytools transcript-assembly -g genome.fasta -i ./clean_data -o ./output
    """

    ta_main = _lazy_import_main()

    args = ['transcript_assembly.py']

    # 必需参数|Required parameters
    args.extend(['-g', genome])
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数|Optional parameters
    if pattern != '*_1.clean.fq.gz':
        args.extend(['-p', pattern])

    if threads != 12:
        args.extend(['-t', str(threads)])

    if sample_timeout != 43200:
        args.extend(['--sample-timeout', str(sample_timeout)])

    if step:
        args.extend(['-s', step])

    if verbose > 0:
        args.extend(['-v'] * verbose)

    if quiet:
        args.append('--quiet')

    if force:
        args.append('--force')

    original_argv = sys.argv
    sys.argv = args

    try:
        ta_main()
    except SystemExit as e:
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
