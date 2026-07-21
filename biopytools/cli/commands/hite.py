"""
HiTE转座子检测与注释命令|HiTE Transposon Detection and Annotation Command
"""

import click
import sys
import os


def _lazy_import_hite_main():
    """
    延迟加载 HiTE 主函数|Lazy load HiTE main function

    Returns:
        callable: HiTE main 函数|HiTE main function

    Raises:
        SystemExit: 导入失败|If import fails
    """
    try:
        from ...hite.main import main as hite_main
        return hite_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    return any(arg in {'-h', '--help'} for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(非帮助模式)|Validate file exists (non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(short_help="HiTE转座子检测与注释|HiTE TE detection and annotation")
@click.option('-i', '--input', 'input',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              help='基因组FASTA文件|Genome FASTA file')
@click.option('-o', '--output-dir', 'output_dir',
              default='./hite_output', show_default=True,
              help='输出目录|Output directory')
@click.option('-t', '--threads', type=int,
              default=12, show_default=True,
              help='线程数|Number of threads')
@click.option('--singularity-path', 'singularity_path',
              default='~/miniforge3/envs/singularity_v.3.8.7/bin/singularity', show_default=True,
              help='Singularity可执行文件|Singularity executable')
@click.option('--sif-file', 'sif_file',
              default='~/software/singularity/hite_3.3.3.sif', show_default=True,
              help='HiTE SIF镜像|HiTE SIF image')
@click.option('--plant/--no-plant', 'plant',
              default=True, show_default=True,
              help='植物基因组|Plant genome')
@click.option('--annotate/--no-annotate', 'annotate',
              default=False, show_default=True,
              help='注释基因组(产出.gff/.out/.tbl)|Annotate genome')
@click.option('--recover/--no-recover', 'recover',
              default=False, show_default=True,
              help='HiTE断点续跑|HiTE recovery mode')
@click.option('--domain/--no-domain', 'domain',
              default=False, show_default=True,
              help='预测TE蛋白结构域|Predict TE domains')
@click.option('--te-type', 'te_type',
              default='all', show_default=True,
              type=click.Choice(['ltr', 'tir', 'helitron', 'non-ltr', 'all']),
              help='TE类型|TE type')
@click.option('--chunk-size', 'chunk_size', type=int,
              default=400, show_default=True,
              help='基因组分块MB|Genome chunk size MB')
@click.option('--miu', type=float,
              default=1.3e-8, show_default=True,
              help='中性突变率|Neutral mutation rate')
@click.option('--min-te-len', 'min_te_len', type=int,
              default=80, show_default=True,
              help='最小TE长度bp|Min TE length bp')
@click.option('--remove-nested/--keep-nested', 'remove_nested',
              default=True, show_default=True,
              help='移除嵌套TE|Remove nested TE')
@click.option('--curated-lib', 'curated_lib',
              default=None,
              help='可信curated TE库(预屏蔽)|Trusted curated TE library')
@click.option('--debug/--no-debug', 'debug',
              default=False, show_default=True,
              help='HiTE debug模式(保留临时文件)|Debug mode')
def hite(input, output_dir, threads, singularity_path, sif_file,
         plant, annotate, recover, domain, te_type, chunk_size, miu,
         min_te_len, remove_nested, curated_lib, debug):
    """
    HiTE转座子检测与注释|HiTE TE detection and annotation

    通过 singularity 直接挂载调用 HiTE,输出到 <output_dir>/01_hite/

    示例|Examples: biopytools hite -i genome.fa -o hite_results -t 12
    """

    hite_main = _lazy_import_hite_main()

    # 构建参数(显式透传所有 default,符合 cli-wrapper-defaults 约定)
    # Build args (explicitly pass all defaults per cli-wrapper-defaults)
    # 注意:main.py 的 bool 参数是 type=int choices=[0,1],必须传 --plant 1 形式
    # Note: main.py bool args are type=int choices=[0,1], must pass --plant 1 form
    args = [
        'hite.py',
        '-i', input,
        '-o', output_dir,
        '-t', str(threads),
        '--singularity-path', singularity_path,
        '--sif-file', sif_file,
        '--plant', '1' if plant else '0',
        '--annotate', '1' if annotate else '0',
        '--recover', '1' if recover else '0',
        '--domain', '1' if domain else '0',
        '--te-type', te_type,
        '--chunk-size', str(chunk_size),
        '--miu', str(miu),
        '--min-te-len', str(min_te_len),
        '--remove-nested', '1' if remove_nested else '0',
        '--debug', '1' if debug else '0',
    ]
    if curated_lib:
        args.extend(['--curated-lib', curated_lib])

    original_argv = sys.argv
    sys.argv = args

    try:
        sys.exit(hite_main())
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
