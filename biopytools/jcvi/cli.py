"""
JCVI CLI组定义|JCVI CLI Group Definition

biopytools jcvi {mcscan,allelic,macro,micro}
"""

import sys
import click

from .mcscan.runner import McscanRunner
from .mcscan.config import McscanConfig
from .allelic.runner import AllelicRunner
from .allelic.config import AllelicConfig


@click.group(
    short_help='JCVI共线性分析工具集|JCVI Synteny Analysis Toolkit',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120),
    invoke_without_command=True,
)
@click.option('--conda-env', default='JCVI_v.1.5.6',
              help='JCVI conda环境名|JCVI conda env name (default: JCVI_v.1.5.6)')
@click.pass_context
def jcvi(ctx, conda_env):
    """JCVI共线性分析工具集 - MCscan/等位基因/宏观共线性/微观共线性|JCVI Synteny Toolkit

    \b
    示例|Examples:
        biopytools jcvi mcscan   -i data -o output
        biopytools jcvi allelic  -i data -o output
        biopytools jcvi macro    -i output --species A,B,C
        biopytools jcvi micro    -i data --pairs A,B --region-a chr1:1-100000 --region-b chr3:1-100000 -o output
    """
    ctx.ensure_object(dict)
    ctx.obj['conda_env'] = conda_env
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@jcvi.command(short_help='MCscan共线性分析|MCscan Collinearity Analysis')
@click.option('-i', '--input', required=True,
              help='[DIR] 输入目录(含*.fa和*.gff)|Input directory with .fa and .gff')
@click.option('-o', '--output', required=True,
              help='[DIR] 输出目录|Output directory')
@click.option('-t', '--threads', type=int, default=24,
              help='[INT] 线程数 (default: 24)')
@click.option('--dbtype', type=click.Choice(['prot', 'nucl']), default='prot',
              help='[STR] 序列类型 (default: prot)')
@click.option('--cscore', type=float, default=0.7,
              help='[FLOAT] C-score阈值 (default: 0.7)')
@click.option('--align-soft', type=click.Choice(['last', 'diamond_blastp']), default='last',
              help='[STR] 比对软件 (default: last)')
@click.option('--pairs', default=None,
              help='[STR] 指定配对 如 A,B A,C|Specific pairs e.g. A,B')
@click.pass_context
def mcscan(ctx, **kwargs):
    """MCscan共线性分析 - 批量两两比较|MCscan Collinearity - Batch pairwise analysis

    示例|Examples: biopytools jcvi mcscan -i data -o output
    """
    config = McscanConfig(
        conda_env=ctx.obj['conda_env'],
        input_dir=kwargs['input'],
        output_dir=kwargs['output'],
        threads=kwargs['threads'],
        dbtype=kwargs['dbtype'],
        cscore=kwargs['cscore'],
        align_soft=kwargs['align_soft'],
        pairs=kwargs['pairs'].split() if kwargs['pairs'] else None,
    )
    runner = McscanRunner(config)
    success = runner.run()
    sys.exit(0 if success else 1)


@jcvi.command(short_help='等位基因鉴定|Allelic Gene Identification')
@click.option('-i', '--input', required=True,
              help='[DIR] 输入目录(含*.fa和*.gff)|Input directory with .fa and .gff')
@click.option('-o', '--output', required=True,
              help='[DIR] 输出目录|Output directory')
@click.option('-t', '--threads', type=int, default=24,
              help='[INT] 线程数 (default: 24)')
@click.option('--dbtype', type=click.Choice(['prot', 'nucl']), default='prot',
              help='[STR] 序列类型 (default: prot)')
@click.option('--cscore', type=float, default=0.7,
              help='[FLOAT] C-score阈值 (default: 0.7)')
@click.option('--align-soft', type=click.Choice(['last', 'diamond_blastp']), default='last',
              help='[STR] 比对软件 (default: last)')
@click.option('--pairs', default=None,
              help='[STR] 指定配对 如 A,B A,C|Specific pairs e.g. A,B')
@click.pass_context
def allelic(ctx, **kwargs):
    """等位基因批量鉴定|Batch Allelic Gene Identification

    示例|Examples: biopytools jcvi allelic -i data -o output
    """
    config = AllelicConfig(
        conda_env=ctx.obj['conda_env'],
        input_dir=kwargs['input'],
        output_dir=kwargs['output'],
        threads=kwargs['threads'],
        dbtype=kwargs['dbtype'],
        cscore=kwargs['cscore'],
        align_soft=kwargs['align_soft'],
        pairs=kwargs['pairs'].split() if kwargs['pairs'] else None,
    )
    runner = AllelicRunner(config)
    success = runner.run()
    sys.exit(0 if success else 1)


@jcvi.command(short_help='宏观共线性可视化|Macrosynteny Karyotype Visualization')
@click.option('-i', '--input', required=True,
              help='[DIR] 输入目录(含*.fa和*.gff)|Input directory with .fa and .gff')
@click.option('-o', '--output', required=True,
              help='[DIR] 输出目录|Output directory')
@click.option('--species', required=True,
              help='[STR] 有序物种列表 A,B,C,... (2=两物种, 3+=多物种)|Ordered species list')
@click.option('-t', '--threads', type=int, default=24,
              help='[INT] 线程数 (default: 24)')
@click.option('--dbtype', type=click.Choice(['prot', 'nucl']), default='prot',
              help='[STR] 序列类型 (default: prot)')
@click.option('--cscore', type=float, default=0.7,
              help='[FLOAT] C-score阈值 (default: 0.7)')
@click.option('--align-soft', type=click.Choice(['last', 'diamond_blastp']), default='last',
              help='[STR] 比对软件 (default: last)')
@click.option('--gff-key', default='ID',
              help='[STR] GFF属性键 (default: ID)')
@click.option('--minspan', type=int, default=30,
              help='[INT] 最小跨度 (default: 30)')
@click.option('--min-chr-genes', type=int, default=20,
              help='[INT] 染色体最小基因数, 少于此的scaffold被过滤 (default: 20)')
@click.option('--figsize', default='',
              help='[STR] 画布大小 如14x10, 多物种默认自动计算|Figure size e.g. 14x10')
@click.option('--shadestyle', type=click.Choice(['line', 'gradient', 'solid']),
              default='line',
              help='[STR] 阴影样式 (default: line)')
@click.option('--chrstyle', type=click.Choice(['rect', 'roundrect']),
              default='rect',
              help='[STR] 染色体样式 (default: rect)')
@click.option('--replot', is_flag=True, default=False,
              help='仅重新绘图, 使用已有的seqids和layout|Re-plot only with existing seqids/layout')
@click.pass_context
def macro(ctx, **kwargs):
    """宏观共线性 - 自动运行mcscan+karyotype可视化(支持多物种)|Macrosynteny - Auto mcscan + karyotype visualization

    示例|Examples: biopytools jcvi macro -i data --species A,B -o output
    """
    from .macro.runner import MacroRunner
    from .macro.config import MacroConfig

    config = MacroConfig(
        conda_env=ctx.obj['conda_env'],
        input_dir=kwargs['input'],
        output_dir=kwargs['output'],
        species=kwargs['species'],
        threads=kwargs['threads'],
        dbtype=kwargs['dbtype'],
        cscore=kwargs['cscore'],
        align_soft=kwargs['align_soft'],
        gff_key=kwargs['gff_key'],
        minspan=kwargs['minspan'],
        min_chr_genes=kwargs['min_chr_genes'],
        figsize=kwargs['figsize'],
        shadestyle=kwargs['shadestyle'],
        chrstyle=kwargs['chrstyle'],
        replot=kwargs['replot'],
    )
    runner = MacroRunner(config)
    success = runner.run()
    sys.exit(0 if success else 1)


@jcvi.command(short_help='微观共线性可视化|Microsynteny Block Visualization')
@click.option('-i', '--input', required=True,
              help='[DIR] 输入目录(含*.fa和*.gff)|Input directory with .fa and .gff')
@click.option('-o', '--output', required=True,
              help='[DIR] 输出目录|Output directory')
@click.option('--pairs', required=True,
              help='[STR] 物种对 A,B|Species pair e.g. A,B')
@click.option('--region-a', required=True,
              help='[STR] 基因组A区间 chr:start-end|Genome A region chr:start-end')
@click.option('--region-b', required=True,
              help='[STR] 基因组B区间 chr:start-end|Genome B region chr:start-end')
@click.option('--genes-a', default=None,
              help='[FILE] 基因组A展示基因文件, 一行一个基因ID|Gene list file for genome A, one per line')
@click.option('--genes-b', default=None,
              help='[FILE] 基因组B展示基因文件, 一行一个基因ID|Gene list file for genome B, one per line')
@click.option('-t', '--threads', type=int, default=12,
              help='[INT] 线程数 (default: 12)')
@click.option('--dbtype', type=click.Choice(['prot', 'nucl']), default='prot',
              help='[STR] 序列类型 (default: prot)')
@click.option('--cscore', type=float, default=0.7,
              help='[FLOAT] C-score阈值 (default: 0.7)')
@click.option('--align-soft', type=click.Choice(['last', 'diamond_blastp']), default='last',
              help='[STR] 比对软件 (default: last)')
@click.option('--glyph-style', type=click.Choice(['arrow', 'box']),
              default='arrow', help='[STR] 基因glyph样式 (default: arrow)')
@click.option('--shadestyle', type=click.Choice(['line', 'gradient', 'solid']),
              default=None, help='[STR] 阴影样式 (default: line)')
@click.pass_context
def micro(ctx, **kwargs):
    """微观共线性 - 自动运行mcscan+区块可视化|Microsynteny - Auto mcscan + block visualization

    示例|Examples: biopytools jcvi micro -i data -o output --pairs A,B --region-a chr1:10000-50000 --region-b chr3:20000-80000
    """
    from .micro.runner import MicroRunner
    from .micro.config import MicroConfig

    config = MicroConfig(
        conda_env=ctx.obj['conda_env'],
        input_dir=kwargs['input'],
        output_dir=kwargs['output'],
        pairs=[kwargs['pairs']],
        region_a=kwargs['region_a'],
        region_b=kwargs['region_b'],
        genes_a=kwargs['genes_a'],
        genes_b=kwargs['genes_b'],
        threads=kwargs['threads'],
        dbtype=kwargs['dbtype'],
        cscore=kwargs['cscore'],
        align_soft=kwargs['align_soft'],
        glyph_style=kwargs['glyph_style'],
        shadestyle=kwargs['shadestyle'] or "",
    )
    runner = MicroRunner(config)
    success = runner.run()
    sys.exit(0 if success else 1)
