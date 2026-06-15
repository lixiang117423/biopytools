"""K-mer工具集命令|K-mer Tools Suite Command"""

import sys
import click


def print_help():
    """打印格式化的帮助信息|Print formatted help"""
    help_text = """
Usage: biopytools kmertools [OPTIONS] COMMAND [ARGS]...

  K-mer工具集 - 建库、查询和分析|K-mer Tools - Build, Query and Analysis

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  build       构建k-mer数据库|Build k-mer database
  compare     比较两个kmer矩阵|Compare two kmer matrix files
  count       K-mer丰度分析|K-mer abundance analysis
  extract     从FASTA提取k-mer|Extract k-mers from FASTA
  gen-fof     生成FOF文件|Generate FOF file
  import-db   导入RocksDB|Import to RocksDB
  intersect   提取目标kmer丰度|Extract target kmer abundance
  kmer2vcf    Kmer丰度转VCF|Kmer abundance to VCF converter
  query       查询k-mer数据库|Query k-mer database
  split-fasta 分割FASTA文件|Split FASTA file
"""
    click.echo(help_text)


@click.command(context_settings=dict(ignore_unknown_options=True, allow_extra_args=True, allow_interspersed_args=False, help_option_names=[]))
@click.pass_context
def kmertools(ctx):
    """K-mer工具集 - 建库、查询和分析|K-mer Tools - Build, Query and Analysis"""
    # 只在没有参数时显示帮助|Show help only when no args
    if not ctx.args:
        print_help()
        return

    # 如果第一个参数是--help/-h且没有子命令，显示帮助|Show help if first arg is --help/-h with no subcommand
    if len(ctx.args) == 1 and ctx.args[0] in ['--help', '-h']:
        print_help()
        return

    try:
        from ...kmertools.main import main as kmertools_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)

    # 执行主程序，传递所有额外参数|Execute main program with all extra arguments
    original_argv = sys.argv
    try:
        # 保留原始参数（包括--help等）|Keep original arguments (including --help, etc.)
        # 跳过'biopytools'和'kmertools'|Skip 'biopytools' and 'kmertools'
        sys.argv = ['kmertools'] + original_argv[2:]
        kmertools_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv


# 导出kmertools命令作为命令入口
__all__ = ['kmertools']
