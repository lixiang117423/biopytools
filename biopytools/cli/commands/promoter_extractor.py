"""
启动子提取器CLI包装器|Promoter Extractor CLI Wrapper
"""

import click
from pathlib import Path
import sys
import os

# 添加biopytools到路径|Add biopytools to path
biopytools_path = Path(__file__).parent.parent.parent
sys.path.insert(0, str(biopytools_path))

from biopytools.promoter_extractor.main import PromoterRunner


@click.command()
@click.option('-g', '--gff', required=True, type=click.Path(exists=True),
              help='输入GFF3文件路径|Input GFF3 file path')
@click.option('--genome', required=True, type=click.Path(exists=True),
              help='输入基因组FASTA文件路径|Input genome FASTA file path', metavar='FASTA')
@click.option('-o', '--output', default='promoters',
              help='输出前缀|Output prefix (default: promoters)')
@click.option('-p', '--promoter-length', default=2000, type=int,
              help='启动子长度（bp）|Promoter length in bp (default: 2000)')
@click.option('--min-length', default=0, type=int,
              help='最小接受长度（bp）|Minimum acceptable length in bp (default: 0)')
@click.option('--gene-list', type=click.Path(exists=True),
              help='基因ID列表文件|Gene ID list file (one gene ID per line)')
@click.option('--no-bed', is_flag=True,
              help='不输出BED格式文件|Do not output BED format file')
@click.option('--no-stats', is_flag=True,
              help='不输出统计文件|Do not output statistics file')
@click.option('-t', '--threads', default=1, type=int,
              help='线程数|Number of threads (default: 1)')
@click.option('--verbose', '-v', count=True,
              help='详细输出模式|Verbose mode (-v: INFO, -vv: DEBUG)')
@click.option('--quiet', is_flag=True,
              help='静默模式(只输出ERROR)|Quiet mode (ERROR only)')
@click.option('--force', '-f', is_flag=True,
              help='强制覆盖已存在文件|Force overwrite existing files')
@click.option('--dry-run', is_flag=True,
              help='模拟运行(不实际执行)|Dry run without execution')
@click.version_option(version='1.0.0')
def promoter_extractor(gff, genome, output, promoter_length, min_length, gene_list,
                       no_bed, no_stats, threads, verbose, quiet, force, dry_run):
    """
    启动子提取器工具|Promoter Extractor Tool

    从GFF3注释文件和基因组FASTA文件中提取基因启动子序列|Extract gene promoter sequences from GFF3 annotation and genome FASTA files

    示例|Examples: biopytools promoter_extractor -g genes.gff3 --genome genome.fa -o promoters
    """
    # 确定日志级别|Determine log level
    if verbose >= 2:
        log_level = "DEBUG"
    elif verbose == 1:
        log_level = "INFO"
    elif quiet:
        log_level = "ERROR"
    else:
        log_level = "INFO"

    try:
        # 创建运行器并运行|Create runner and run
        runner = PromoterRunner(
            gff_file=gff,
            genome_file=genome,
            output_prefix=output,
            promoter_length=promoter_length,
            min_length=min_length,
            gene_list=gene_list,
            output_bed=not no_bed,
            output_stats=not no_stats,
            threads=threads,
            keep_intermediate=False,
            log_level=log_level,
            quiet=quiet,
            verbose=verbose,
            force=force,
            dry_run=dry_run
        )

        if dry_run:
            click.echo("模拟运行模式-不会实际执行命令|DRY RUN mode - commands will not be executed")

        # 执行提取|Run extraction
        runner.run()

    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        raise click.ClickException(str(e))


if __name__ == '__main__':
    promoter_extractor()
