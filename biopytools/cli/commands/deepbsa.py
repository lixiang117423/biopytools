"""
DeepBSA批量分析命令|DeepBSA Batch Analysis Command
支持子命令模式：batch, run, merge|Supports subcommands: batch, run, merge
"""

import click
import sys


def _lazy_import_deepbsa_cli():
    """延迟加载deepbsa CLI主函数|Lazy load deepbsa CLI main function"""
    try:
        from ...deepbsa.cli import main as deepbsa_cli
        return deepbsa_cli
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


@click.command(
    short_help='DeepBSA批量分析工具（支持子命令）|DeepBSA batch analysis tool (with subcommands)',
    context_settings=dict(help_option_names=[], ignore_unknown_options=True, max_content_width=120),
    name='deepbsa'
)
@click.argument('args', nargs=-1, type=click.UNPROCESSED)
def deepbsa(args):
    """
    DeepBSA批量分析工具（支持子命令）|DeepBSA Batch Analysis Tool (with subcommands)

    子命令|Subcommands:
      batch   生成批量处理命令（推荐）|Generate batch processing commands (Recommended)
      run     运行单个DeepBSA分析方法|Run single DeepBSA analysis method
      merge   合并DeepBSA运行结果|Merge DeepBSA results
      vcf2csv VCF转CSV（为DeepBSA准备输入数据）|Convert VCF to CSV for DeepBSA

    推荐工作流|Recommended workflow:
      1. biopytools deepbsa batch -i input.vcf -o batch_jobs/
      2. [投递并运行 batch_jobs/run.sh]
      3. bash batch_jobs/all/merge_results.sh

    示例|Examples:
      biopytools deepbsa batch -i input.vcf -o batch_jobs/
      biopytools deepbsa run -i input.vcf -o results/ -m DL,K
      biopytools deepbsa merge -i results/ -o merged/

    可用方法|Available methods: DL, K, ED4, SNP, SmoothG, SmoothLOD, Ridit

    运行 "biopytools deepbsa batch -h" 查看batch子命令的详细帮助
    Run "biopytools deepbsa batch -h" for batch subcommand detailed help
    """
    # 延迟加载|Lazy loading
    deepbsa_cli_main = _lazy_import_deepbsa_cli()

    # 构建参数列表|Build argument list
    sys.argv = ['deepbsa'] + list(args)

    try:
        deepbsa_cli_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)

