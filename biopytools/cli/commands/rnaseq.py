"""
RNA-seq分析命令|RNA-seq Analysis Command
"""

import click
import sys
import os


def _lazy_import_rnaseq_main():
    """延迟加载rnaseq主函数|Lazy load rnaseq main function"""
    try:
        from ...rnaseq.main import main as rnaseq_main
        return rnaseq_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='RNA-seq分析流程|RNA-seq analysis pipeline',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组FASTA文件路径|Genome FASTA file path')
@click.option('--gtf', '-f',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因注释GTF文件路径|Gene annotation GTF file path')
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='输入FASTQ文件目录或样本信息文件|Input FASTQ file directory or sample information file')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--pattern', '-p',
              type=str,
              help='FASTQ文件命名模式|Fastq file naming pattern (e.g., "*.R1.fastq.gz" or "*_1.fq.gz"), * represents sample name')
@click.option('--remove', '-r',
              default='no',
              type=click.Choice(['yes', 'y', 'no', 'n']),
              show_default=True,
              help='处理后删除BAM文件|Remove BAM files after processing')
@click.option('--verbose', '-v',
              count=True,
              help='增加输出详细程度|Increase output verbosity (-v, -vv, -vvv)')
@click.option('--quiet',
              is_flag=True,
              help='静默模式|Quiet mode, only output errors')
@click.option('--log-file',
              type=str,
              help='日志文件路径|Log file path')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
              default='INFO',
              show_default=True,
              help='日志级别|Log level')
@click.option('--dry-run',
              is_flag=True,
              help='试运行模式|Dry run mode, no actual execution')
@click.option('--threads', '-t',
              type=int,
              default=8,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--sample-timeout',
              type=int,
              default=21600,
              show_default=True,
              help='单个样本处理超时时间（秒）|Sample processing timeout in seconds')
def rnaseq(genome, gtf, input, output, pattern, remove, verbose, quiet,
           log_file, log_level, dry_run, threads, sample_timeout):
    """
    RNA-seq分析流程：HISAT2 + StringTie|RNA-seq Analysis Pipeline: HISAT2 + StringTie

    完整的RNA-seq数据分析流程，包括HISAT2比对和StringTie定量，计算FPKM和TPM值|Complete RNA-seq data analysis pipeline, including HISAT2 alignment and StringTie quantification, calculating FPKM and TPM values

    示例|Examples:
    biopytools rnaseq -g genome.fa -f genes.gtf -i /data/fastq/ -o rnaseq_results
    """

    # 延迟加载|Lazy load
    rnaseq_main = _lazy_import_rnaseq_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['rnaseq.py']

    # 必需参数|Required parameters
    args.extend(['-g', genome])
    args.extend(['-f', gtf])
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if pattern is not None:
        args.extend(['-p', pattern])

    if remove != 'no':
        args.extend(['-r', remove])

    if threads != 8:
        args.extend(['-t', str(threads)])

    if sample_timeout != 21600:
        args.extend(['--sample-timeout', str(sample_timeout)])

    # 日志选项|Logging options
    if verbose > 0:
        args.extend(['-v'] * verbose)

    if quiet:
        args.append('--quiet')

    if log_file is not None:
        args.extend(['--log-file', log_file])

    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    # 高级选项|Advanced options
    if dry_run:
        args.append('--dry-run')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        rnaseq_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
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
