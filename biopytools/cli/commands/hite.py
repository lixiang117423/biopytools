"""
HiTE转座子检测与注释命令|HiTE Transposon Detection and Annotation Command
"""

import click
import sys
import os


def _lazy_import_hite_main():
    """
    延迟加载HiTE主函数|Lazy load HiTE main function

    Returns:
        callable: HiTE主函数|HiTE main function

    Raises:
        SystemExit: 如果导入失败|If import fails
    """
    try:
        from ...hite.hite_analyzer import main as hite_main
        return hite_main
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
    short_help='HiTE转座子检测与注释|HiTE transposon detection and annotation',
    context_settings=dict(
        help_option_names=['-h', '--help'],
        max_content_width=120
    )
)
@click.option(
    '--input', '-i',
    required=True,
    type=click.Path(exists=True),
    callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
    help='基因组FASTA文件|Genome FASTA file path'
)
@click.option(
    '--output-dir', '-o',
    default='./hite_output',
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
    '--plant/--no-plant',
    default=True,
    show_default=True,
    help='是否为植物基因组|Whether the genome is from plant'
)
@click.option(
    '--annotate',
    is_flag=True,
    help='是否注释基因组|Whether to annotate the genome with detected TEs'
)
@click.option(
    '--recover',
    is_flag=True,
    help='是否启用断点续跑|Whether to enable recovery mode'
)
@click.option(
    '--domain',
    is_flag=True,
    help='是否预测TE蛋白结构域|Whether to predict protein domains in TEs'
)
@click.option(
    '--te-type',
    default='all',
    type=click.Choice(['ltr', 'tir', 'helitron', 'non-ltr', 'all']),
    help='检测的TE类型|TE type to detect'
)
@click.option(
    '--chunk-size',
    default=400,
    show_default=True,
    type=int,
    help='基因组分块大小(MB)|Genome chunk size in MB'
)
@click.option(
    '--miu',
    default=1.3e-8,
    show_default=True,
    type=float,
    help='中性突变率(per bp per year)|Neutral mutation rate'
)
@click.option(
    '--min-te-len',
    default=80,
    show_default=True,
    type=int,
    help='最小TE长度(bp)|Minimum TE length in bp'
)
def hite(input, output_dir, threads, singularity_cmd, sif_file,
         plant, annotate, recover, domain, te_type, chunk_size, miu, min_te_len):
    """
    HiTE转座子检测与注释工具|HiTE Transposon Detection and Annotation Tool

    使用HiTE通过Singularity容器进行基因组转座子检测与注释
    Use HiTE through Singularity container for genome transposon detection and annotation

    \b
    示例|Examples: biopytools hite -i genome.fa -o hite_results -t 12
    """

    # 延迟加载|Lazy loading
    hite_main = _lazy_import_hite_main()

    # 构建参数列表|Build argument list
    kwargs = {
        'genome': input,  # Note: parameter name is 'input' but config expects 'genome'
        'output_dir': output_dir,
        'threads': threads,
        'singularity_cmd': singularity_cmd,
        'sif_file': sif_file,
        'plant': plant,
        'annotate': annotate,
        'recover': recover,
        'domain': domain,
        'te_type': te_type,
        'chunk_size': chunk_size,
        'miu': miu,
        'min_te_len': min_te_len
    }

    # 执行主程序|Execute main program
    try:
        analyzer = hite_main(**kwargs)
        analyzer.run_analysis()
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
