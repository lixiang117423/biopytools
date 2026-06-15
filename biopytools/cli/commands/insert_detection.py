"""插入检测CLI包装器|Insert detection CLI wrapper"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from biopytools.insert_detection.main import main as insert_detection_main
        return insert_detection_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在|Validate path exists"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(short_help="插入序列位点检测|Insert sequence insertion site detection")
@click.option('-i', '--genome',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='参考基因组FASTA文件|Reference genome FASTA file')
@click.option('--insert',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='插入序列FASTA文件|Insert sequence FASTA file')
@click.option('--fastq-dir',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='FASTQ文件目录|FASTQ files directory')
@click.option('-o', '--output-dir',
              default='./insert_detection_output',
              help='输出目录|Output directory (default: ./insert_detection_output)')
@click.option('-t', '--threads',
              default=12,
              show_default=True,
              help='线程数|Threads')
@click.option('--min-clip',
              type=int,
              default=20,
              show_default=True,
              help='最小soft-clip长度|Minimum soft-clip length')
@click.option('--min-support',
              type=int,
              default=5,
              show_default=True,
              help='最小支持reads数|Minimum supporting reads')
@click.option('--score-threshold',
              type=int,
              default=1000,
              show_default=True,
              help='得分阈值|Score threshold')
@click.option('--bowtie2-path',
              default='bowtie2',
              show_default=True,
              help='Bowtie2可执行文件路径|Bowtie2 executable path')
@click.option('--samtools-path',
              default='samtools',
              show_default=True,
              help='samtools可执行文件路径|samtools executable path')
@click.option('--read1-suffix',
              default='_1.clean.fq.gz',
              show_default=True,
              help='R1文件后缀（包含扩展名）|Read 1 file suffix with extension')
@click.option('--read2-suffix',
              default='_2.clean.fq.gz',
              show_default=True,
              help='R2文件后缀（包含扩展名）|Read 2 file suffix with extension')
@click.option('--force',
              is_flag=True,
              help='强制重新运行所有步骤|Force rerun all steps')
@click.option('--verbose',
              is_flag=True,
              help='显示详细日志|Verbose logging')
@click.option('--quiet',
              is_flag=True,
              help='仅显示错误|Errors only')
def insert_detection(genome, insert, fastq_dir, output_dir, threads,
                     min_clip, min_support, score_threshold,
                     bowtie2_path, samtools_path, read1_suffix, read2_suffix,
                     force, verbose, quiet):
    """
    插入序列位点检测流程|Insert sequence insertion site detection pipeline

    检测插入序列在基因组中的插入位置，自动识别配对的FASTQ文件|Detect insertion sites of insert sequence in genome, auto-identify paired FASTQ files

    默认使用biopytools fastp输出格式（_1.clean.fq.gz / _2.clean.fq.gz），支持自定义R1/R2后缀|Default uses biopytools fastp output format, supports custom R1/R2 suffixes

    示例|Examples: biopytools insert-detection -i genome.fa --insert tdna.fa --fastq-dir fastq_output/ -o output/
    """
    insert_detection_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['insert_detection.py']
    args.extend(['-i', genome])
    args.extend(['--insert', insert])
    args.extend(['--fastq-dir', fastq_dir])
    args.extend(['-o', output_dir])
    args.extend(['-t', str(threads)])
    args.extend(['--min-clip', str(min_clip)])
    args.extend(['--min-support', str(min_support)])
    args.extend(['--score-threshold', str(score_threshold)])
    args.extend(['--bowtie2-path', bowtie2_path])
    args.extend(['--samtools-path', samtools_path])
    args.extend(['--read1-suffix', read1_suffix])
    args.extend(['--read2-suffix', read2_suffix])
    if force:
        args.append('--force')
    if verbose:
        args.append('--verbose')
    if quiet:
        args.append('--quiet')

    original_argv = sys.argv
    sys.argv = args

    try:
        insert_detection_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
