"""
转录组验证注释|Transcriptome Validation for Genome Annotation Command
"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...rnaseq_val.main import main as rnaseq_val_main
        return rnaseq_val_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path exists (only in non-help mode)"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='转录组验证注释|Transcriptome validation for genome annotation',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-g', '--genome', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='基因组FASTA文件|Genome FASTA file')
@click.option('-a', '--annotation', required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='参考注释GTF文件|Reference annotation GTF file')
@click.option('-o', '--output', required=True,
              help='输出目录|Output directory')
@click.option('--sr-dir', default=None,
              help='二代clean reads目录|SR clean reads directory')
@click.option('--lr-dir', default=None,
              help='三代clean reads目录|LR clean reads directory')
@click.option('--lr-platform', default='pacbio_hifi',
              type=click.Choice(['pacbio_hifi', 'ont']),
              help='三代测序平台|LR platform')
@click.option('-t', '--threads', default=16,
              help='线程数|Threads')
@click.option('--strandness', default='RF',
              type=click.Choice(['RF', 'FR', 'unstranded']),
              help='文库链特异性|Library strandness')
@click.option('--steps', default='all',
              help='执行步骤(逗号分隔或all)|Steps (comma-separated or all)')
@click.option('--force', is_flag=True,
              help='强制重新运行|Force re-run')
def rnaseq_val(genome, annotation, output, sr_dir, lr_dir, lr_platform,
               threads, strandness, steps, force):
    """
    转录组验证注释工具|Transcriptome validation for genome annotation

    示例|Examples: biopytools rnaseq-val -g genome.fa -a anno.gtf --sr-dir ./sr/ -o out
    """

    rnaseq_val_main = _lazy_import_main()

    # 构建参数列表|Build argument list
    args = ['rnaseq_val.py']

    args.extend(['-g', genome])
    args.extend(['-a', annotation])
    args.extend(['-o', output])

    if sr_dir:
        args.extend(['--sr-dir', sr_dir])
    if lr_dir:
        args.extend(['--lr-dir', lr_dir])
    if lr_platform != 'pacbio_hifi':
        args.extend(['--lr-platform', lr_platform])
    if threads != 16:
        args.extend(['-t', str(threads)])
    if strandness != 'RF':
        args.extend(['--strandness', strandness])
    if steps != 'all':
        args.extend(['--steps', steps])
    if force:
        args.append('--force')

    original_argv = sys.argv
    sys.argv = args

    try:
        rnaseq_val_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
