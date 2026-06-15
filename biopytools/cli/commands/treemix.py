"""TreeMix CLI适配器|TreeMix CLI Adapter

将 biopytools treemix 注册到LazyGroup, 采用jcvi透传模式
Routes biopytools treemix to internal click.Group via sys.argv passthrough
"""

import sys
import click


def print_help():
    help_text = """\
Usage: biopytools treemix COMMAND [OPTIONS]

  TreeMix群体历史与基因流分析|TreeMix Population History & Gene Flow

Options:
  -h, --help  Show this message and exit.

Commands:
  prepare  输入准备 - VCF转TreeMix格式|Input Preparation
  run      运行TreeMix|Run TreeMix
  all      完整流程 - 输入准备+运行|Full Pipeline
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
def treemix(ctx):
    """TreeMix群体历史与基因流分析|TreeMix Population History & Gene Flow"""
    if not ctx.args:
        print_help()
        return

    if ctx.args[0] in ['--help', '-h']:
        print_help()
        return

    try:
        from ...treemix.cli import treemix as treemix_group
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)

    original_argv = sys.argv
    try:
        sys.argv = ['treemix'] + ctx.args
        treemix_group(standalone_mode=False)
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv


__all__ = ['treemix']
