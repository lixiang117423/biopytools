"""
Kraken2宏基因组分类命令|Kraken2 Metagenomic Classification Command
"""

import click
import sys
import os


def _lazy_import_kraken2_main():
    """延迟加载kraken2主函数|Lazy load kraken2 main function"""
    try:
        from ...kraken2.main import main as kraken2_main
        return kraken2_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_dir_exists(dir_path):
    """验证目录存在(仅在非帮助模式)|Validate directory exists (only in non-help mode)"""
    if not _is_help_request() and dir_path and not os.path.isdir(dir_path):
        raise click.BadParameter(f"目录不存在|Directory does not exist: {dir_path}")
    return dir_path


@click.command(
    short_help='Kraken2宏基因组分类工具|Kraken2 metagenomic classification tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input-dir', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_dir_exists(value) if value else None,
              help='输入FASTQ目录|Input FASTQ directory')
@click.option('--db-path', '-d',
              required=True,
              callback=lambda ctx, param, value: _validate_dir_exists(value) if value else None,
              help='Kraken2数据库路径|Kraken2 database path')
@click.option('--output-dir', '-o',
              default='./kraken2_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              default=12,
              show_default=True,
              type=int,
              help='线程数|Number of threads')
@click.option('--read-len',
              default=150,
              show_default=True,
              type=int,
              help='读长(用于Bracken)|Read length (for Bracken)')
@click.option('--confidence',
              default=0.0,
              show_default=True,
              type=float,
              help='Kraken2置信度阈值|Kraken2 confidence score threshold')
@click.option('--bracken-level',
              default='S',
              show_default=True,
              type=click.Choice(['D', 'P', 'C', 'O', 'F', 'G', 'S', 'S1']),
              help='Bracken分类级别|Bracken taxonomic level')
@click.option('--bracken-threshold',
              default=10,
              show_default=True,
              type=int,
              help='Bracken最小读数阈值|Bracken minimum read threshold')
@click.option('--no-bracken',
              is_flag=True,
              help='跳过Bracken分析|Skip Bracken analysis')
@click.option('--r1-suffix',
              default='_1.clean.fq.gz',
              show_default=True,
              help='R1文件后缀|R1 file suffix')
@click.option('--r2-suffix',
              default='_2.clean.fq.gz',
              show_default=True,
              help='R2文件后缀|R2 file suffix')
def kraken2(input_dir, db_path, output_dir, threads, read_len, confidence,
            bracken_level, bracken_threshold, no_bracken, r1_suffix, r2_suffix):
    """
    Kraken2宏基因组分类工具|Kraken2 Metagenomic Classification Tool

    自动识别配对FASTQ并使用Kraken2进行物种分类|Auto-detect paired FASTQ and run Kraken2 species classification

    示例|Examples: biopytools kraken2 -i ./fastq/ -d ~/database/kraken2_db -o ./kraken2_output
    """

    kraken2_main = _lazy_import_kraken2_main()

    args = ['kraken2.py']
    args.extend(['-i', input_dir])
    args.extend(['-d', db_path])

    if output_dir != './kraken2_output':
        args.extend(['-o', output_dir])
    if threads != 12:
        args.extend(['-t', str(threads)])
    if read_len != 150:
        args.extend(['--read-len', str(read_len)])
    if confidence != 0.0:
        args.extend(['--confidence', str(confidence)])
    if bracken_level != 'S':
        args.extend(['--bracken-level', bracken_level])
    if bracken_threshold != 10:
        args.extend(['--bracken-threshold', str(bracken_threshold)])
    if r1_suffix != '_1.clean.fq.gz':
        args.extend(['--r1-suffix', r1_suffix])
    if r2_suffix != '_2.clean.fq.gz':
        args.extend(['--r2-suffix', r2_suffix])
    if no_bracken:
        args.append('--no-bracken')

    original_argv = sys.argv
    sys.argv = args

    try:
        kraken2_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
