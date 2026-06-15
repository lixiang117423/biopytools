"""
EDTA转座子注释命令|EDTA TE Annotation Command
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


@click.command(short_help='EDTA转座子注释|EDTA TE Annotation',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('-i', '--genome', required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组FASTA文件|Genome FASTA file')
@click.option('--species', default='others', show_default=True,
              type=click.Choice(['Rice', 'Maize', 'others']),
              help='物种类型|Species type')
@click.option('--step', default='all', show_default=True,
              type=click.Choice(['all', 'filter', 'final', 'anno']),
              help='运行步骤|Step to run')
@click.option('--overwrite', default=0, show_default=True,
              type=click.IntRange(0, 1),
              help='覆盖已有结果|Overwrite existing results')
@click.option('--cds',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='CDS序列文件|CDS sequences file')
@click.option('--curatedlib',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='筛选TE库|Curated TE library')
@click.option('--rmlib',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='RepeatModeler库|RepeatModeler library')
@click.option('--sensitive', default=0, show_default=True,
              type=click.IntRange(0, 1),
              help='使用RepeatModeler|Use RepeatModeler')
@click.option('--anno', default=0, show_default=True,
              type=click.IntRange(0, 1),
              help='执行全基因组注释|Perform whole-genome annotation')
@click.option('--rmout',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='RepeatMasker输出文件|RepeatMasker output file')
@click.option('--maxdiv', default=40, show_default=True,
              type=int,
              help='最大分歧度|Maximum divergence')
@click.option('--evaluate', default=0, show_default=True,
              type=click.IntRange(0, 1),
              help='评估注释一致性|Evaluate annotation consistency')
@click.option('--exclude',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='排除区域BED文件|Exclude regions BED file')
@click.option('--force', default=0, show_default=True,
              type=click.IntRange(0, 1),
              help='强制使用水稻TE|Force to use rice TEs')
@click.option('--u', default=1.3e-8, show_default=True,
              type=float,
              help='中性突变率|Neutral mutation rate')
@click.option('-t', '--threads', default=12, show_default=True,
              type=int,
              help='线程数|Number of threads')
@click.option('--debug', default=0, show_default=True,
              type=click.IntRange(0, 1),
              help='保留中间文件|Retain intermediate files')
@click.option('-o', '--output-dir', default='./edta_output', show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--edta-path',
              type=click.Path(),
              help='EDTA安装路径|EDTA installation path')
def edta(genome, species, step, overwrite, cds, curatedlib, rmlib,
         sensitive, anno, rmout, maxdiv, evaluate, exclude, force, u,
         threads, debug, output_dir, edta_path):
    """
    EDTA转座子注释工具|EDTA TE Annotation Tool

    单基因组转座子鉴定、分类和注释|Single-genome transposon identification, classification, and annotation

    示例|Examples: biopytools edta -i plant.fa --anno 1 -t 24
    """
    # 延迟加载|Lazy loading
    edta_main = _lazy_import_edta_main()

    # 构建参数列表|Build argument list
    args = ['edta.py', 'edta']

    # 必需参数|Required parameters
    args.extend(['--genome', genome])

    # 可选参数|Optional parameters
    if species != 'others':
        args.extend(['--species', species])

    if step != 'all':
        args.extend(['--step', step])

    if overwrite != 0:
        args.extend(['--overwrite', str(overwrite)])

    if cds:
        args.extend(['--cds', cds])

    if curatedlib:
        args.extend(['--curatedlib', curatedlib])

    if rmlib:
        args.extend(['--rmlib', rmlib])

    if sensitive != 0:
        args.extend(['--sensitive', str(sensitive)])

    if anno != 0:
        args.extend(['--anno', str(anno)])

    if rmout:
        args.extend(['--rmout', rmout])

    if maxdiv != 40:
        args.extend(['--maxdiv', str(maxdiv)])

    if evaluate != 0:
        args.extend(['--evaluate', str(evaluate)])

    if exclude:
        args.extend(['--exclude', exclude])

    if force != 0:
        args.extend(['--force', str(force)])

    if u != 1.3e-8:
        args.extend(['--u', str(u)])

    if threads != 12:
        args.extend(['--threads', str(threads)])

    if debug != 0:
        args.extend(['--debug', str(debug)])

    if output_dir != './edta_output':
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

