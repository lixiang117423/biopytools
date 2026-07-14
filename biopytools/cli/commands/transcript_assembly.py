"""
转录本组装命令|Transcript Assembly Command
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
    if path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='转录本组装(FASTQ或BAM→GFF3)|Transcript assembly (FASTQ or BAM → GFF3)',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--genome', '-g',
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='基因组FASTA(FASTQ模式或--transcripts时必需)|Genome FASTA (required for FASTQ mode or --transcripts)')
@click.option('--input', '-i',
              type=click.Path(exists=True),
              help='输入FASTQ文件目录(与-b互斥)|Input FASTQ dir, mutually exclusive with -b')
@click.option('--bam', '-b',
              multiple=True,
              type=click.Path(exists=True),
              help='输入BAM文件(可多次-b,与-i互斥)|Input BAM file(s), repeatable, mutually exclusive with -i')
@click.option('--guide-gff',
              type=click.Path(exists=True),
              help='参考注释GTF/GFF3(-G guided)|Reference annotation for guided assembly')
@click.option('--read-type',
              type=click.Choice(['auto', 'short', 'long']),
              default='auto', show_default=True,
              help='读长类型|Read type')
@click.option('--transcripts',
              is_flag=True,
              help='额外输出transcripts.fa(需-g)|Also output cDNA (needs -g)')
@click.option('--predict-cds',
              is_flag=True,
              help='TransDecoder预测CDS(需-g,输出gene/mRNA/CDS)|TransDecoder CDS prediction (needs -g)')
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
              type=click.Choice(['1', '2', '3', '4', '5', '6', '7']),
              help='运行指定步骤|Run only specified step')
@click.option('--verbose', '-v',
              count=True,
              help='增加输出详细程度|Increase output verbosity')
@click.option('--quiet',
              is_flag=True,
              help='静默模式，仅输出错误信息|Quiet mode, only output errors')
@click.option('--force',
              is_flag=True,
              help='强制重新处理已完成的步骤|Force re-process completed steps')
def transcript_assembly(output, genome, input, bam, guide_gff, read_type, transcripts, predict_cds,
                         pattern, threads, sample_timeout, step, verbose, quiet, force):
    """
    转录本组装流程:HISAT2/StringTie,支持FASTQ或BAM直入、短/长读、GFF3输出|Transcript assembly: HISAT2/StringTie, FASTQ or BAM input, short/long reads, GFF3 output

    示例|Examples: biopytools transcript-assembly -b sample.bam -o ./out
    """

    ta_main = _lazy_import_main()

    args = ['transcript_assembly.py']

    # 必需参数|Required parameters
    args.extend(['-o', output])

    # 输入参数(条件透传)|Input parameters (conditional passthrough)
    if genome:
        args.extend(['-g', genome])
    if input:
        args.extend(['-i', input])
    if bam:
        args.append('-b')
        args.extend(bam)
    if guide_gff:
        args.extend(['--guide-gff', guide_gff])
    if read_type != 'auto':
        args.extend(['--read-type', read_type])
    if transcripts:
        args.append('--transcripts')
    if predict_cds:
        args.append('--predict-cds')

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
