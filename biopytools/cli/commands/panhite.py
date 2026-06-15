"""
panHiTE群体基因组TE分析命令|panHiTE Pan-genome TE Analysis Command
"""

import click
import sys
import os


def _lazy_import_panhite_main():
    """
    延迟加载panHiTE主函数|Lazy load panHiTE main function

    Returns:
        callable: panHiTE主函数|panHiTE main function

    Raises:
        SystemExit: 如果导入失败|If import fails
    """
    try:
        from ...hite.panhite_analyzer import main as panhite_main
        return panhite_main
    except ImportError as e:
        click.echo(
            f"导入错误|Import Error: {e}",
            err=True
        )
        sys.exit(1)


def _is_help_request():
    """
    检查是否为帮助请求|Check if this is a help request

    Returns:
        bool: 如果是帮助请求返回True|Returns True if help request
    """
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """
    验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)

    Args:
        file_path: 文件路径|File path

    Returns:
        str: 文件路径|File path

    Raises:
        click.BadParameter: 如果文件不存在|If file does not exist
    """
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(
            f"文件不存在|File does not exist: {file_path}"
        )
    return file_path


@click.command(
    short_help='panHiTE群体基因组TE分析|panHiTE pan-genome TE analysis',
    context_settings=dict(
        help_option_names=['-h', '--help'],
        max_content_width=120
    )
)
@click.option(
    '--pan-genomes-dir', '-p',
    required=True,
    type=click.Path(exists=True),
    callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
    help='泛基因组目录(包含所有基因组文件)|Pan-genomes directory containing all genome files'
)
@click.option(
    '--input', '-i',
    required=True,
    type=click.Path(exists=True),
    callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
    help='基因组列表文件|Genome list file path'
)
@click.option(
    '--genes-dir',
    type=click.Path(exists=True),
    callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
    help='基因注释文件目录(GFF/GTF)|Genes annotation files directory (GFF/GTF)'
)
@click.option(
    '--rna-dir',
    type=click.Path(exists=True),
    callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
    help='RNA-seq数据目录|RNA-seq data directory'
)
@click.option(
    '--output-dir', '-o',
    default='./panhite_output',
    show_default=True,
    type=click.Path(),
    help='输出目录|Output directory path'
)
@click.option(
    '--threads', '-t',
    default=12,
    show_default=True,
    type=int,
    help='线程数|Number of threads'
)
@click.option(
    '--singularity-cmd',
    default='~/miniforge3/envs/singularity_v.3.8.7/bin/singularity',
    show_default=True,
    help='Singularity命令路径|Singularity executable path'
)
@click.option(
    '--sif-file',
    default='~/software/singularity/hite_3.3.3.sif',
    show_default=True,
    help='HiTE SIF镜像路径|HiTE SIF image path'
)
@click.option(
    '--te-type',
    default='all',
    type=click.Choice(['ltr', 'tir', 'helitron', 'non-ltr', 'all']),
    help='检测的TE类型|TE type to detect'
)
@click.option(
    '--miu',
    default=1.3e-8,
    show_default=True,
    type=float,
    help='中性突变率(per bp per year)|Neutral mutation rate'
)
@click.option(
    '--skip-analyze',
    is_flag=True,
    help='跳过分析，仅生成panTE库|Skip analysis, only generate panTE library'
)
@click.option(
    '--recover',
    is_flag=True,
    help='是否启用断点续跑|Whether to enable recovery mode'
)
@click.option(
    '--debug',
    is_flag=True,
    help='是否开启调试模式|Whether to enable debug mode'
)
def panhite(pan_genomes_dir, input, genes_dir, rna_dir, output_dir,
            threads, singularity_cmd, sif_file, te_type, miu,
            skip_analyze, recover, debug):
    """
    panHiTE群体基因组转座子分析工具|panHiTE Pan-genome Transposon Analysis Tool

    使用panHiTE通过Singularity容器进行泛基因组转座子检测与比较分析
    Use panHiTE through Singularity container for pan-genome transposon detection
    and comparative analysis

    \b
    基本用法|Basic usage:
        biopytools panhite -p pan_genomes/ -i genome_list.txt -o results

    \b
    包含基因表达分析|With gene expression analysis:
        biopytools panhite -p pan_genomes/ -i genome_list.txt \\
            --genes-dir annotations/ --rna-dir rnaseq/ -o results

    \b
    仅生成panTE库|Generate panTE library only:
        biopytools panhite -p pan_genomes/ -i genome_list.txt --skip-analyze

    \b
    基因组列表文件格式|Genome list file format:
        # 基本格式(仅TE检测)|Basic format (TE detection only)
        genome1.fa
        genome2.fa

        # 包含基因注释|With gene annotation
        genome1.fa    annotation1.gff
        genome2.fa    annotation2.gff

        # 包含RNA-seq(用于TIDELs分析)|With RNA-seq (for TIDELs analysis)
        genome1.fa    annotation1.gff    1    rna1_1.fq.gz    rna1_2.fq.gz
        genome2.fa    annotation2.gff    1    rna2_1.fq.gz    rna2_2.fq.gz

    \b
    参考|Reference:
        Hu, K. et al. HiTE: accurate transposable element identification and
        pan-genome analysis. bioRxiv preprint. https://doi.org/10.1101/2025.02.15.638472v1
    """

    # 延迟加载|Lazy loading
    panhite_main = _lazy_import_panhite_main()

    # 构建参数列表|Build argument list
    kwargs = {
        'pan_genomes_dir': pan_genomes_dir,
        'genome_list': input,  # Note: parameter name is 'input' but config expects 'genome_list'
        'output_dir': output_dir,
        'threads': threads,
        'singularity_cmd': singularity_cmd,
        'sif_file': sif_file,
        'te_type': te_type,
        'miu': miu,
        'skip_analyze': skip_analyze,
        'recover': recover,
        'debug': debug
    }

    # 可选参数|Optional parameters
    if genes_dir:
        kwargs['genes_dir'] = genes_dir
    if rna_dir:
        kwargs['rna_dir'] = rna_dir

    # 执行主程序|Execute main program
    try:
        analyzer = panhite_main(**kwargs)
        analyzer.run_analysis()
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
