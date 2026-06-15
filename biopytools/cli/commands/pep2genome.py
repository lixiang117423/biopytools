"""
蛋白质到基因组比对命令|Protein to Genome Alignment Command
"""

import click
import sys
import os


def _lazy_import_pep2genome_main():
    """延迟加载pep2genome主函数|Lazy load pep2genome main function"""
    try:
        from ...pep2genome.main import main as pep2genome_main
        return pep2genome_main
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
    short_help='蛋白质到基因组比对工具|Protein to genome alignment tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              help='基因组FASTA文件|Genome FASTA file')
@click.option('--protein', '-p',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              help='蛋白质FASTA文件|Protein FASTA file')
@click.option('--output-dir', '-o',
              default='./pep2genome_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              default=12,
              show_default=True,
              type=int,
              help='线程数|Number of threads')
@click.option('--miniprot-path',
              default='miniprot',
              show_default=True,
              help='Miniprot工具路径|Miniprot tool path')
@click.option('--no-gff3',
              is_flag=True,
              help='不导出GFF3格式|Do not export GFF3 format')
@click.option('--no-bed',
              is_flag=True,
              help='不导出BED格式|Do not export BED format')
@click.option('--no-statistics',
              is_flag=True,
              help='不生成统计报告|Do not generate statistics report')
def pep2genome(genome, protein, output_dir, threads, miniprot_path, no_gff3, no_bed, no_statistics):
    """
    蛋白质到基因组比对工具|Protein to Genome Alignment Tool

    使用Miniprot进行蛋白质与基因组的比对，支持PAF格式解析、统计分析和GFF3/BED格式导出
    Align proteins to genome using Miniprot, with PAF parsing, statistics analysis, and GFF3/BED export

    示例|Example: biopytools pep2genome --genome genome.fa --protein protein.fa -o results
    """

    # 延迟加载|Lazy loading
    pep2genome_main = _lazy_import_pep2genome_main()

    # 构建参数列表|Build argument list
    args = ['pep2genome.py']

    # 必需参数|Required parameters
    args.extend(['--genome', genome])
    args.extend(['--protein', protein])
    args.extend(['-o', output_dir])

    # 可选参数|Optional parameters
    if threads != 12:
        args.extend(['-t', str(threads)])

    if miniprot_path != 'miniprot':
        args.extend(['--miniprot-path', miniprot_path])

    # 输出选项|Output options
    if no_gff3:
        args.append('--no-gff3')

    if no_bed:
        args.append('--no-bed')

    if no_statistics:
        args.append('--no-statistics')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        pep2genome_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
