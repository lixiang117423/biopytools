"""
Merqury QV计算命令|Merqury QV Calculation Command
"""

import click
import sys
import os


def _lazy_import_assembly_qv_main():
    """延迟加载assembly_qv主函数|Lazy load assembly_qv main function"""
    try:
        from ...assembly_qv.main import main as assembly_qv_main
        return assembly_qv_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path exists (only in non-help mode)"""
    if not _is_help_request():
        if not os.path.exists(path):
            raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='装配质量QV计算|Assembly Quality QV Calculation',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='FASTQ文件或目录|FASTQ file or directory')
@click.option('-g', '--genome',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='基因组FASTA文件|Genome FASTA file')
@click.option('-o', '--output-dir',
              default='./assembly_qv_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('-k', '--kmer-size',
              type=int,
              default=None,
              help='K-mer大小|K-mer size (default: auto-select)')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--conda-env',
              default="~/miniforge3/envs/merqury_v.1.3/bin/",
              show_default=True,
              help='Conda环境路径|Conda environment path')
@click.option('--data-type',
              type=click.Choice(['auto', 'illumina', 'hifi']),
              default='auto',
              show_default=True,
              help='数据类型|Data type')
def assembly_qv(input, genome, output_dir, kmer_size, threads, conda_env, data_type):
    """
    装配质量QV值计算工具|Assembly Quality QV Value Calculation Tool

    使用k-mer光谱分析计算基因组装配的QV值|Calculate QV value for genome assembly using k-mer spectrum analysis

    示例|Examples: biopytools assembly-qv -i fastq_dir/ -g genome.fa
    """

    # 延迟加载|Lazy loading
    assembly_qv_main = _lazy_import_assembly_qv_main()

    # 构建参数列表|Build argument list
    args = ['assembly_qv.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-g', genome])

    # 可选参数|Optional parameters
    if output_dir != './assembly_qv_output':
        args.extend(['-o', output_dir])

    if kmer_size is not None:
        args.extend(['-k', str(kmer_size)])

    if threads != 24:
        args.extend(['-t', str(threads)])

    if conda_env != "~/miniforge3/envs/merqury_v.1.3/bin/":
        args.extend(['--conda-env', conda_env])

    if data_type != 'auto':
        args.extend(['--data-type', data_type])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        assembly_qv_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
