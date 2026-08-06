"""PSVCP 泛基因组构建命令|PSVCP pangenome construction command"""

import os
import sys

import click


def _lazy_import_runner():
    """延迟加载 PSVCPRunner|Lazy load PSVCPRunner"""
    try:
        from ...psvcp.main import PSVCPRunner
        return PSVCPRunner
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    return any(arg in {'-h', '--help'} for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(非帮助模式)|Validate path exists (non-help mode)"""
    if not _is_help_request() and path and not os.path.exists(os.path.expanduser(path)):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='PSVCP 线性泛基因组构建(MUMmer+Assemblytics)|PSVCP linear pangenome construction',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--genome-dir',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              type=click.Path(),
              help='genome 目录(含 {name}.fa + {name}.gff/.gff3)|genome dir with {name}.fa + {name}.gff/.gff3')
@click.option('-l', '--genome-list',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              type=click.Path(),
              help='genome_list 文本(行1=ref,其余=query)|genome_list (line1=ref, rest=queries)')
@click.option('-o', '--output-dir',
              default='~/psvcp_out',
              show_default=True,
              type=click.Path(),
              help='输出目录|output directory')
@click.option('-t', '--threads',
              default=12,
              show_default=True,
              type=int,
              help='线程数|threads')
@click.option('--force',
              is_flag=True,
              help='忽略断点续传,强制重跑|ignore checkpoint, force rerun')
@click.option('--log-file',
              type=click.Path(),
              help='日志文件路径(默认 output_dir/psvcp.log)|log file path')
@click.option('-v', '--verbose',
              is_flag=True,
              help='详细输出|verbose output')
def psvcp(genome_dir, genome_list, output_dir, threads, force, log_file, verbose):
    """PSVCP 线性泛基因组构建|PSVCP linear pangenome construction

    基于多个组装基因组,用 MUMmer+Assemblytics 检测 PAV 并依次并入参考,构建线性泛基因组。
    |Build a linear pan-genome by incorporating PAVs (MUMmer+Assemblytics) from each query genome.

    示例|Examples: biopytools psvcp -i genome_gff_dir/ -l genome_list.txt -o pangenome_out/
    """

    try:
        PSVCPRunner = _lazy_import_runner()
        runner = PSVCPRunner(
            genome_dir=str(genome_dir),
            genome_list=str(genome_list),
            output_dir=str(output_dir),
            threads=threads,
            force=force,
            log_file=str(log_file) if log_file else None,
            verbose=verbose,
        )
        success = runner.run()
        if success:
            click.echo("PSVCP 泛基因组构建完成|PSVCP pangenome construction completed!")
        else:
            click.echo("PSVCP 泛基因组构建失败|PSVCP pangenome construction failed!", err=True)
            sys.exit(1)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"发生错误|Error occurred: {e}", err=True)
        sys.exit(1)
