"""
基因组组装质量评估命令|Genome Assembly Quality Control Command

biopytools assembly-qc命令的Click包装器
Click wrapper for biopytools assembly-qc command
"""

import click
import sys
import os


def _lazy_import_assembly_qc_main():
    """延迟加载assembly_qc主函数|Lazy load assembly_qc main function"""
    try:
        from ...assembly_qc.main import main as assembly_qc_main
        return assembly_qc_main
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
    name='assembly-qc',
    short_help='基因组组装质量评估|Genome Assembly Quality Control',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组FASTA文件|Genome FASTA file')
@click.option('--lineage', '-l',
              required=True,
              help='BUSCO数据集名称（如embryophyta_odb10）或完整路径|BUSCO dataset name (e.g., embryophyta_odb10) or full path')
@click.option('--output-dir', '-o',
              default='./assembly_qc_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--sample-name', '-s',
              default='genome_sample',
              show_default=True,
              help='样品名称|Sample name')
@click.option('--skip-busco',
              is_flag=True,
              help='跳过BUSCO评估|Skip BUSCO evaluation')
@click.option('--skip-lai',
              is_flag=True,
              help='跳过LAI评估|Skip LAI evaluation')
@click.option('--lai-full-mode',
              is_flag=True,
              help='LAI完整模式（不使用-qq，运行blastn计算，用于种间比较）|LAI full mode (no -qq, run blastn for interspecies comparison)')
@click.option('--enable-qv',
              is_flag=True,
              help='启用QV评估（默认启用）|Enable QV evaluation (default: enabled)')
@click.option('--qv-kmer-size',
              type=int,
              help='k-mer大小（None表示自动选择）|K-mer size (None for auto)')
@click.option('--enable-mapping',
              is_flag=True,
              help='启用NGS Mapping评估（默认启用）|Enable NGS mapping evaluation (default: enabled)')
@click.option('--enable-long-read-mapping',
              is_flag=True,
              help='启用三代数据Mapping评估（默认启用）|Enable long-read mapping evaluation (default: enabled)')
@click.option('--mapping-pattern',
              default='_1.clean.fq.gz',
              show_default=True,
              help='FASTQ文件匹配模式|FASTQ file pattern')
@click.option('--ngs-reads',
              type=click.Path(exists=True),
              help='NGS reads目录（用于QV和mapping）|NGS reads directory (for QV and mapping)')
@click.option('--long-reads',
              type=click.Path(exists=True),
              help='Long-reads目录（用于QV和mapping）|Long-reads directory (for QV and mapping)')
@click.option('--long-read-type',
              type=click.Choice(['ont', 'pacbio', 'hifi']),
              default='hifi',
              show_default=True,
              help='Long-read数据类型|Long-read data type')
@click.option('--no-html',
              is_flag=True,
              help='不生成HTML报告|Do not generate HTML report')
@click.option('--no-table',
              is_flag=True,
              help='不生成表格|Do not generate table')
@click.option('--table-format',
              type=click.Choice(['tsv', 'xlsx', 'both']),
              default='both',
              show_default=True,
              help='表格格式|Table format')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数（自动分配给各子模块）|Threads (automatically distributed to sub-modules)')
def assembly_qc(genome, lineage, output_dir, sample_name,
                skip_busco, skip_lai, lai_full_mode,
                enable_qv, qv_kmer_size,
                enable_mapping, enable_long_read_mapping, mapping_pattern,
                ngs_reads, long_reads, long_read_type,
                no_html, no_table, table_format, threads):
    """
    基因组组装质量评估工具|Genome Assembly Quality Control Tool

    提供基因组组装质量的综合评估，包括BUSCO完整性、LAI指数、QV值和Mapping统计
    Provides comprehensive genome assembly quality evaluation, including BUSCO completeness, LAI index, QV value, and mapping statistics

    示例|Examples: biopytools assembly-qc --genome genome.fa --lineage embryophyta_odb10 --ngs-reads ./illumina_reads --long-reads ./hifi_reads --long-read-type hifi -o qc_results
    """

    # 延迟加载|Lazy loading
    assembly_qc_main_func = _lazy_import_assembly_qc_main()

    # 构建参数列表|Build argument list
    args = [
        'assembly_qc.py',
        '--genome', genome,
        '--lineage', lineage,
        '--output-dir', output_dir,
        '--sample-name', sample_name,
    ]

    # 核心评估|Core evaluation
    if skip_busco:
        args.append('--skip-busco')
    if skip_lai:
        args.append('--skip-lai')
    if lai_full_mode:
        args.append('--lai-full-mode')

    # QV评估|QV evaluation
    if enable_qv:
        args.append('--enable-qv')
        if qv_kmer_size:
            args.extend(['--qv-kmer-size', str(qv_kmer_size)])

    # Mapping评估|Mapping evaluation
    if enable_mapping:
        args.append('--enable-mapping')
        args.extend(['--mapping-pattern', mapping_pattern])

    # 三代数据Mapping评估|Long-read mapping evaluation
    if enable_long_read_mapping:
        args.append('--enable-long-read-mapping')

    # Reads数据|Reads data
    if ngs_reads:
        args.extend(['--ngs-reads', ngs_reads])
    if long_reads:
        args.extend(['--long-reads', long_reads])
        args.extend(['--long-read-type', long_read_type])

    # 报告参数|Report parameters
    if no_html:
        args.append('--no-html')
    if no_table:
        args.append('--no-table')
    args.extend(['--table-format', table_format])

    # 线程参数|Thread parameter
    args.extend(['--threads', str(threads)])

    # 运行主程序|Run main program
    sys.argv = args
    exit_code = assembly_qc_main_func()
    sys.exit(exit_code)
