"""
LAI模块CLI包装器|LAI Module CLI Wrapper
"""

import click
import sys


def _lazy_import_lai_main():
    """延迟加载LAI主函数|Lazy load LAI main function"""
    try:
        from ...lai.main import LAICalculator
        return LAICalculator
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


@click.command()
@click.option('-i', '--input', required=True, type=click.Path(exists=True),
              help='基因组FASTA文件|Genome FASTA file')
@click.option('-o', '--output', required=True,
              help='输出目录|Output directory')
@click.option('-t', '--threads', default=12, show_default=True,
              type=int,
              help='线程数|Number of threads')
@click.option('-m', '--mode', default='full', show_default=True,
              type=click.Choice(['full', 'harvest', 'retrieve', 'calculate']),
              help='运行模式: full(完整流程), harvest(仅候选识别), retrieve(仅筛选), calculate(仅LAI计算)|Run mode')
@click.option('--skip-completed/--no-skip-completed', default=True, show_default=True,
              is_flag=True,
              help='跳过已完成的步骤|Skip completed steps')
@click.option('--conda-harvest',
              default='~/miniforge3/envs/ltr_harvest_parallel_v.1.2',
              show_default=True,
              help='LTR_harvest conda环境路径|LTR_harvest conda environment path')
@click.option('--conda-finder',
              default='~/miniforge3/envs/ltr_finder_parallel_v.1.3',
              show_default=True,
              help='LTR_finder conda环境路径|LTR_finder conda environment path')
@click.option('--conda-retriever',
              default='~/miniforge3/envs/ltr_retriever_v.3.0.1',
              show_default=True,
              help='LTR_retriever conda环境路径|LTR_retriever conda environment path')
def lai(input, output, threads, mode, skip_completed, conda_harvest, conda_finder, conda_retriever):
    """
    LAI计算工具|LTR Assembly Index Calculator

    用于评估基因组组装质量的长末端重复逆转座子组装指数计算工具
    Tool for calculating LTR Assembly Index to evaluate genome assembly quality

    示例|Examples: biopytools lai -i genome.fa -o output_dir
    """
    try:
        # 延迟加载|Lazy load
        LAICalculator = _lazy_import_lai_main()

        # 创建计算器并运行|Create calculator and run
        calculator = LAICalculator(
            genome=input,
            output_dir=output,
            threads=threads,
            mode=mode,
            skip_completed=skip_completed,
            conda_env_ltr_harvest=conda_harvest,
            conda_env_ltr_finder=conda_finder,
            conda_env_ltr_retriever=conda_retriever
        )

        calculator.run()

    except Exception as e:
        click.echo(click.style(f"错误|Error: {str(e)}", fg='red'), err=True)
        sys.exit(1)
