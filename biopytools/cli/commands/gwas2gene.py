"""
GWAS候选基因筛选命令|GWAS Candidate Gene Finder Command
"""

import click
import sys


def _lazy_import_gwas2gene_main():
    """延迟加载gwas2gene主函数|Lazy load gwas2gene main function"""
    try:
        from ...gwas2gene.main import main as gwas2gene_main
        return gwas2gene_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


@click.command(
    short_help='GWAS候选基因筛选工具|GWAS candidate gene finder tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--gwas', '-g',
              required=True,
              type=click.Path(exists=True),
              help='GWAS结果文件路径|GWAS result file path')
@click.option('--pval-col', '-p',
              required=True,
              type=str,
              help='P值所在列名或列号（1-based）|P-value column name or index (1-based)')
@click.option('--threshold', '-t',
              default=1e-5,
              type=float,
              show_default=True,
              help='P值阈值|P-value threshold')
@click.option('--window', '-w',
              default=200000,
              type=int,
              show_default=True,
              help='上下游窗口大小|Window size upstream/downstream in bp')
@click.option('--gff',
              required=True,
              type=click.Path(exists=True),
              help='GFF3注释文件路径|GFF3 annotation file path')
@click.option('--output', '-o',
              required=True,
              type=str,
              help='输出文件路径|Output file path')
@click.option('--func', '-f',
              type=click.Path(exists=True),
              help='功能注释文件路径|Function annotation file path (optional)')
def gwas2gene(gwas, pval_col, threshold, window, gff, output, func):
    """
    GWAS候选基因筛选工具|GWAS Candidate Gene Finder Tool

    根据GWAS结果筛选显著SNP，在指定窗口内从GFF文件中提取候选基因|Extract candidate genes near significant GWAS SNPs from GFF annotations

    示例|Examples: biopytools gwas2gene -g gwas.txt -p Pvalue -a annotation.gff3 -o candidates.tsv
    """

    # 延迟加载|Lazy load
    gwas2gene_main_func = _lazy_import_gwas2gene_main()

    # 构建参数列表|Build argument list
    args = ['gwas2gene.py']

    # 必需参数|Required parameters
    args.extend(['--gwas', gwas])
    args.extend(['--pval-col', pval_col])
    args.extend(['--gff', gff])
    args.extend(['--output', output])

    # 可选参数（只在非默认时添加）|Optional parameters (add only when non-default)
    if threshold != 1e-5:
        args.extend(['--threshold', str(threshold)])

    if window != 200000:
        args.extend(['--window', str(window)])

    if func:
        args.extend(['--func', func])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        gwas2gene_main_func()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
