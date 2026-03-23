"""
互作转录组分析命令|Dual RNA-seq Analysis Command
"""

import click
import sys
import os


def _lazy_import_dual_rnaseq_main():
    """延迟加载dual_rnaseq主函数|Lazy load dual_rnaseq main function"""
    try:
        from ...dual_rnaseq.main import main as dual_rnaseq_main
        return dual_rnaseq_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
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


@click.command(
    short_help='互作转录组分析工具|Dual RNA-seq analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--species1-name',
              required=True,
              help='物种1名称|Species 1 name (e.g., host)')
@click.option('--species1-genome',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              help='物种1基因组FASTA文件|Species 1 genome FASTA file')
@click.option('--species1-gtf',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              help='物种1GTF注释文件|Species 1 GTF annotation file')
@click.option('--species2-name',
              required=True,
              help='物种2名称|Species 2 name (e.g., pathogen)')
@click.option('--species2-genome',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              help='物种2基因组FASTA文件|Species 2 genome FASTA file')
@click.option('--species2-gtf',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              help='物种2GTF注释文件|Species 2 GTF annotation file')
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTQ文件目录或样本信息文件|Input FASTQ file directory or sample information file')
@click.option('--output-dir', '-o',
              default='./dual_rnaseq_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--pattern', '-p',
              default='*_1.clean.fq.gz',
              show_default=True,
              help='FASTQ文件命名模式|FASTQ file naming pattern')
@click.option('--threads', '-t',
              default=12,
              show_default=True,
              type=int,
              help='线程数|Number of threads')
@click.option('--min-mapq',
              default=20,
              show_default=True,
              type=int,
              help='最小mapping quality值|Minimum mapping quality value')
@click.option('--no-unique-only',
              is_flag=True,
              help='不禁用非唯一比对|Do not disable non-unique mappings')
@click.option('--no-extract-fastq',
              is_flag=True,
              help='不提取FASTQ文件|Do not extract FASTQ files from BAM (default: extract)')
def dual_rnaseq(species1_name, species1_genome, species1_gtf,
               species2_name, species2_genome, species2_gtf,
               input, output_dir, pattern, threads, min_mapq, no_unique_only, no_extract_fastq):
    """
    互作转录组分析工具|Dual RNA-seq Analysis Tool

    双物种转录组同时分析，包括物种分类、定量和表达矩阵生成|Simultaneous analysis of dual-species transcriptome, including species classification, quantification and expression matrix generation

    示例|Example: biopytools dual-rnaseq --species1-name host --species1-genome host.fa --species1-gtf host.gtf --species2-name pathogen --species2-genome pathogen.fa --species2-gtf pathogen.gtf -i ./fastq_data -o results
    """

    # 延迟加载|Lazy loading
    dual_rnaseq_main = _lazy_import_dual_rnaseq_main()

    # 构建参数列表|Build argument list
    args = ['dual_rnaseq.py']

    # 物种1参数|Species 1 parameters
    args.extend(['--species1-name', species1_name])
    args.extend(['--species1-genome', species1_genome])
    args.extend(['--species1-gtf', species1_gtf])

    # 物种2参数|Species 2 parameters
    args.extend(['--species2-name', species2_name])
    args.extend(['--species2-genome', species2_genome])
    args.extend(['--species2-gtf', species2_gtf])

    # 通用参数|Common parameters
    args.extend(['-i', input])

    if output_dir != './dual_rnaseq_output':
        args.extend(['-o', output_dir])

    if pattern != '*_1.clean.fq.gz':
        args.extend(['-p', pattern])

    if threads != 12:
        args.extend(['-t', str(threads)])

    if min_mapq != 20:
        args.extend(['--min-mapq', str(min_mapq)])

    if no_unique_only:
        args.append('--no-unique-only')

    if no_extract_fastq:
        args.append('--no-extract-fastq')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        dual_rnaseq_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
