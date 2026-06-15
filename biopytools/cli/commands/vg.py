"""VG变异图分析命令|VG Variation Graph Analysis Command"""

import sys
import click


def print_help():
    """打印格式化的帮助信息|Print formatted help"""
    help_text = """
Usage: biopytools vg [OPTIONS] COMMAND [ARGS]...

  VG变异图分析工具|VG Variation Graph Analysis Tool

Options:
  -h, --help     Show this message and exit.

Commands:
  construct      从VCF和参考基因组构建变异图|Build variation graph from VCF and reference
  deconstruct    从变异图生成VCF|Export VCF from variation graph
  giraffe        快速序列比对|Fast read alignment
  index          创建图索引|Create graph indexes
"""
    click.echo(help_text)


@click.command(context_settings=dict(ignore_unknown_options=True,
                                     allow_extra_args=True,
                                     allow_interspersed_args=False,
                                     help_option_names=[]))
@click.pass_context
def vg(ctx):
    """VG变异图分析工具|VG Variation Graph Analysis Tool"""
    # 只在没有参数时显示帮助|Show help only when no args
    if not ctx.args:
        print_help()
        return

    # 如果第一个参数是--help/-h且没有子命令，显示帮助|Show help if first arg is --help/-h with no subcommand
    if len(ctx.args) == 1 and ctx.args[0] in ['--help', '-h']:
        print_help()
        return

    try:
        from ...vg.main import main as vg_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)

    # 执行主程序，传递所有额外参数|Execute main program with all extra arguments
    original_argv = sys.argv
    try:
        # 保留原始参数（包括--help等）|Keep original arguments (including --help, etc.)
        # 跳过'biopytools'和'vg'|Skip 'biopytools' and 'vg'
        sys.argv = ['vg'] + original_argv[2:]
        vg_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv


# 导出vg命令作为命令入口
__all__ = ['vg']
