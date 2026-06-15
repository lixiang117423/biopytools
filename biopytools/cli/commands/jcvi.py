"""
JCVI CLI适配器|JCVI CLI Adapter

将 biopytools jcvi 注册到LazyGroup, 采用kmertools透传模式
Routes biopytools jcvi to internal click.Group via sys.argv passthrough
"""

import sys
import click


def print_help():
    help_text = """
Usage: biopytools jcvi COMMAND [OPTIONS]

  JCVI共线性分析工具集|JCVI Synteny Analysis Toolkit

Options:
  --conda-env TEXT  JCVI conda环境名 (default: JCVI_v.1.5.6)
  -h, --help        Show this message and exit.

Commands:
  mcscan    MCscan共线性分析|MCscan Collinearity Analysis
  allelic   等位基因批量鉴定|Allelic Gene Batch Identification
  macro     宏观共线性可视化|Macrosynteny Karyotype Visualization
  micro     微观共线性可视化|Microsynteny Block Visualization
"""
    click.echo(help_text)


@click.command(
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
        allow_interspersed_args=False,
        help_option_names=[],
    )
)
@click.pass_context
def jcvi(ctx):
    """JCVI共线性分析工具集|JCVI Synteny Analysis Toolkit"""
    if not ctx.args:
        print_help()
        return

    if ctx.args[0] in ['--help', '-h']:
        print_help()
        return

    try:
        from ...jcvi.cli import jcvi as jcvi_group
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)

    original_argv = sys.argv
    try:
        sys.argv = ['jcvi'] + ctx.args
        jcvi_group(standalone_mode=False)
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv


__all__ = ['jcvi']
