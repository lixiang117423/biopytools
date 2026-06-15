"""
PanEDTA泛基因组转座子注释命令|PanEDTA Pan-genome TE Annotation Command
"""

import click
import sys
import os


def _lazy_import_edta_main():
    """延迟加载EDTA主函数|Lazy load EDTA main function"""
    try:
        from ...edta.main import main as edta_main
        return edta_main
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


@click.command(short_help='PanEDTA泛基因组转座子注释|PanEDTA Pan-genome TE Annotation',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('-i', '--genome-list', required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组列表文件|Genome list file')
@click.option('-c', '--cds',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='CDS序列文件|CDS sequences file')
@click.option('-l', '--curatedlib',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='筛选TE库|Curated TE library')
@click.option('-f', '--fl-copy', default=3, show_default=True,
              type=int,
              help='全长拷贝数阈值|Full-length copy number cutoff')
@click.option('-a', '--anno', default=1, show_default=True,
              type=click.IntRange(0, 1),
              help='执行全基因组注释|Perform whole-genome annotation')
@click.option('--overwrite', default=0, show_default=True,
              type=click.IntRange(0, 1),
              help='覆盖已有结果|Overwrite existing results')
@click.option('-t', '--threads', default=12, show_default=True,
              type=int,
              help='线程数|Number of threads')
@click.option('-o', '--output-dir', default='./panedta_output', show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--edta-path',
              type=click.Path(),
              help='EDTA安装路径|EDTA installation path')
def panedta(genome_list, cds, curatedlib, fl_copy, anno, overwrite,
           threads, output_dir, edta_path):
    """
    PanEDTA泛基因组转座子注释工具|PanEDTA Pan-genome TE Annotation Tool

    多个基因组的转座子联合分析和注释|Joint transposon analysis and annotation for multiple genomes

    示例|Examples: biopytools panedta -i genomes.txt -c cds.fa -t 24
    """
    # 延迟加载|Lazy loading
    edta_main = _lazy_import_edta_main()

    # 构建参数列表|Build argument list
    args = ['edta.py', 'panedta']

    # 必需参数|Required parameters
    args.extend(['--genome-list', genome_list])

    # 可选参数|Optional parameters
    if cds:
        args.extend(['-c', cds])

    if curatedlib:
        args.extend(['-l', curatedlib])

    if fl_copy != 3:
        args.extend(['--fl-copy', str(fl_copy)])

    if anno != 1:
        args.extend(['--anno', str(anno)])

    if overwrite != 0:
        args.extend(['--overwrite', str(overwrite)])

    if threads != 12:
        args.extend(['--threads', str(threads)])

    if output_dir != './panedta_output':
        args.extend(['--output-dir', output_dir])

    if edta_path:
        args.extend(['--edta-path', edta_path])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        edta_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断操作|Analysis interrupted by user", err=True)
        sys.exit(130)
    except Exception as e:
        click.echo(f"分析执行失败|Analysis execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
