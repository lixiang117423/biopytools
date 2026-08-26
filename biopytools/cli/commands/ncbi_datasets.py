"""NCBI datasets 批量下载命令|NCBI datasets batch download command"""

import click
import sys
import os


def _lazy_import_ncbi_main():
    """延迟加载ncbi_datasets主函数|Lazy load ncbi_datasets main function"""
    try:
        from ...ncbi_datasets.main import main as ncbi_main
        return ncbi_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


@click.command(short_help='NCBI taxon基因组批量下载|NCBI taxon genome batch download',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--taxon', '-t',
              required=True,
              type=int,
              help='NCBI taxon 编号|NCBI taxon ID')
@click.option('--output-dir', '-o',
              default='./output',
              help='输出目录|Output directory')
@click.option('--assembly-source',
              type=click.Choice(['refseq', 'genbank']),
              help='只下载 RefSeq 或 GenBank 来源|Only RefSeq or GenBank assemblies')
@click.option('--assembly-level',
              help='组装级别过滤(逗号分隔)|Assembly level filter (comma-separated): '
                   'complete, chromosome, scaffold, contig')
@click.option('--reference',
              is_flag=True,
              help='只下载参考基因组|Reference genomes only')
@click.option('--annotated',
              is_flag=True,
              help='只下载有注释的基因组|Annotated genomes only')
@click.option('--include-gff3',
              is_flag=True,
              help='额外下载 GFF3 基因注释|Also download GFF3 gene annotation')
@click.option('--include-protein',
              is_flag=True,
              help='额外下载蛋白序列|Also download protein sequences')
@click.option('--include-cds',
              is_flag=True,
              help='额外下载 CDS 序列|Also download CDS sequences')
@click.option('--include-seq-report',
              is_flag=True,
              help='额外下载 seq-report 汇总|Also download seq-report summary')
@click.option('--dry-run',
              is_flag=True,
              help='只查询 assembly 清单,不下载|Query manifest only, no download')
@click.option('--organize/--no-organize',
              default=True,
              show_default=True,
              help='下载后整理到 02_organized|Organize into 02_organized after download')
@click.option('--datasets-path',
              help='datasets 工具路径(默认走环境变量/配置/~bin)|datasets tool path '
                   '(default: env var / config / ~/bin)')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO',
              help='日志级别|Log level')
def ncbi_datasets(taxon, output_dir, assembly_source, assembly_level, reference,
                  annotated, include_gff3, include_protein, include_cds,
                  include_seq_report, dry_run, organize, datasets_path, log_level):
    """NCBI taxon 基因组批量下载|NCBI taxon genome batch download

    输入 taxon 编号,用官方 datasets CLI 下载该 taxon 下所有基因组
    |Input a taxon ID and download all genomes under it via the official datasets CLI

    示例|Examples: biopytools ncbi-datasets -t 67593 -o output/
    """

    ncbi_main = _lazy_import_ncbi_main()

    # 构建参数列表|Build argument list
    args = ['ncbi_datasets.py']
    args.extend(['-t', str(taxon)])
    if output_dir != './output':
        args.extend(['-o', output_dir])
    if assembly_source:
        args.extend(['--assembly-source', assembly_source])
    if assembly_level:
        args.extend(['--assembly-level', assembly_level])
    if reference:
        args.append('--reference')
    if annotated:
        args.append('--annotated')
    if include_gff3:
        args.append('--include-gff3')
    if include_protein:
        args.append('--include-protein')
    if include_cds:
        args.append('--include-cds')
    if include_seq_report:
        args.append('--include-seq-report')
    if dry_run:
        args.append('--dry-run')
    if not organize:
        args.append('--no-organize')
    if datasets_path:
        args.extend(['--datasets-path', datasets_path])
    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    original_argv = sys.argv
    sys.argv = args

    try:
        ncbi_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
