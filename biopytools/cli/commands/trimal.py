"""trimal 多序列比对修剪命令|trimal MSA Trimming Command"""

import click
import sys
import os


def _lazy_import_trimal_runner():
    """延迟加载 TrimalRunner|Lazy load TrimalRunner"""
    try:
        from ...trimal.main import TrimalRunner
        return TrimalRunner
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


@click.command(
    short_help='多序列比对自动修剪(trimAl)|MSA automated trimming with trimAl',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='输入比对文件|Input alignment file')
@click.option('--output-dir', '-o',
              default='./trimal_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--method',
              default='automated1',
              show_default=True,
              type=click.Choice(['automated1', 'gappyout', 'strict', 'strictplus', 'gt', 'cons']),
              help='修剪方法|Trimming method (automated1 适合 ML 建树|suited for ML trees)')
@click.option('--gt-threshold',
              default=0.9,
              show_default=True,
              type=float,
              help='gap 阈值[0,1],仅 method=gt|gap threshold [0,1], only method=gt')
@click.option('--cons-threshold',
              default=80,
              show_default=True,
              type=int,
              help='保守度阈值[0,100],仅 method=cons|conservation [0,100], only method=cons')
@click.option('--format',
              default='keep',
              show_default=True,
              type=click.Choice(['keep', 'fasta', 'phylip', 'phylip_paml', 'clustal', 'nexus', 'nbrf', 'mega']),
              help='输出格式(keep=沿用输入)|Output format (keep=input format)')
@click.option('--colnumbering',
              is_flag=True,
              help='输出新旧列号映射|Output old/new column mapping')
@click.option('--backtrans',
              type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
              help='CDS 文件,AA→NT 反向翻译|CDS file for AA→NT backtranslation')
@click.option('--complementary',
              is_flag=True,
              help='输出被修剪列的互补比对|Output complementary alignment of trimmed columns')
@click.option('--keep-header',
              is_flag=True,
              help='保留完整 FASTA 头|Keep full FASTA headers')
@click.option('--sample-name',
              help='输出文件名前缀(默认输入 basename)|Output filename prefix (default: input basename)')
@click.option('--log-file',
              type=click.Path(),
              help='日志文件路径|Log file path')
@click.option('--verbose', '-v',
              is_flag=True,
              help='详细输出|Verbose output')
def trimal(input, output_dir, method, gt_threshold, cons_threshold, format,
           colnumbering, backtrans, complementary, keep_header, sample_name,
           log_file, verbose):
    """多序列比对自动修剪|MSA automated trimming with trimAl

    示例|Examples: biopytools trimal -i aln.fasta -o out/ --method automated1
    """

    try:
        TrimalRunner = _lazy_import_trimal_runner()

        runner = TrimalRunner(
            input_file=str(input),
            output_dir=str(output_dir),
            method=method,
            gt_threshold=gt_threshold,
            cons_threshold=cons_threshold,
            output_format=format,
            colnumbering=colnumbering,
            backtrans_file=str(backtrans) if backtrans else None,
            complementary=complementary,
            keep_header=keep_header,
            sample_name=sample_name,
            log_file=str(log_file) if log_file else None,
            verbose=verbose
        )

        success = runner.run_extraction()

        if success:
            click.echo("trimal 修剪完成|trimal trimming completed successfully!")
        else:
            click.echo("trimal 修剪失败|trimal trimming failed!", err=True)
            sys.exit(1)

    except KeyboardInterrupt:
        click.echo("\n用户中断|Interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"发生错误|Error occurred: {e}", err=True)
        sys.exit(1)
