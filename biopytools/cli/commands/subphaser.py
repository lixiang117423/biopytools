"""SubPhaser亚基因组分离工具|SubPhaser subgenome phasing command"""

import click
import os
import sys


def _lazy_import_subphaser_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...subphase.main import main as subphaser_main
        return subphaser_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    return any(arg in {'-h', '--help'} for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在|Validate path exists"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='异源多倍体亚基因组分离|Phase subgenomes of allopolyploids',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)

# 必需参数|Required parameters
@click.option('-i', '--genomes', required=True, multiple=True,
              callback=lambda ctx, param, value: [_validate_path_exists(p) for p in value] if value else None,
              help='基因组FASTA文件(可多个)|Genome FASTA files')

@click.option('--nsg', required=True, type=int,
              help='亚基因组数量(>=2)|Number of subgenomes (>=2)')

# 输入模式|Input mode
@click.option('-c', '--sg-cfgs', multiple=True,
              callback=lambda ctx, param, value: [_validate_path_exists(p) for p in value] if value else None,
              help='亚基因组配置文件(可选，不提供则自动模式)|Subgenome config files (optional, auto mode if omitted)')

@click.option('--parental-genomes', multiple=True,
              callback=lambda ctx, param, value: [_validate_path_exists(p) for p in value] if value else None,
              help='父本基因组(验证模式，需2个)|Parental genomes for validation (requires 2)')

# 输出|Output
@click.option('-o', '--output-dir', default='./subphaser_output',
              show_default=True,
              help='输出目录|Output directory')

@click.option('--prefix', default=None,
              help='输出前缀|Output prefix')

# 资源|Resources
@click.option('-t', '--threads', type=int, default=24,
              show_default=True,
              help='线程数|Number of threads')

# 染色体过滤|Chromosome filtering
@click.option('--min-chrom-size', type=int, default=1_000_000,
              show_default=True,
              help='最小染色体长度(bp)，过滤小contigs|Min chromosome size (bp), filter small contigs')

# K-mer参数|K-mer parameters
@click.option('-k', '--kmer-size', type=int, default=15,
              show_default=True,
              help='K-mer大小|K-mer size')

@click.option('-f', '--min-fold', type=float, default=2.0,
              show_default=True,
              help='最小倍数差异|Minimum fold difference')

@click.option('-q', '--min-freq', type=int, default=200,
              show_default=True,
              help='最小k-mer频率|Minimum k-mer frequency')

# 聚类|Cluster
@click.option('--max-pval', type=float, default=0.05,
              show_default=True,
              help='最大P值|Maximum P-value')

@click.option('--replicates', type=int, default=1000,
              show_default=True,
              help='Bootstrap重复次数|Bootstrap replicates')

@click.option('--test-method', default='ttest_ind',
              type=click.Choice(['ttest_ind', 'kruskal', 'wilcoxon', 'mannwhitneyu']),
              show_default=True,
              help='统计检验方法|Statistical test method')

# 步骤控制|Step control
@click.option('--disable-ltr', is_flag=True,
              help='禁用LTR分析(大基因组耗时较长)|Disable LTR analysis')

@click.option('--disable-circos', is_flag=True,
              help='禁用Circos图|Disable Circos plot')

@click.option('--disable-blocks', is_flag=True,
              help='禁用同源区块|Disable homologous blocks')

@click.option('--just-core', is_flag=True,
              help='仅运行核心phasing|Only run core phasing')

# LTR参数|LTR parameters
@click.option('--ltr-detectors', multiple=True,
              type=click.Choice(['ltr_finder', 'ltr_harvest']),
              help='LTR检测工具|LTR detection tools')

@click.option('--mu', type=float, default=13e-9,
              show_default=True,
              help='替换率/年|Substitution rate per year')

# Circos参数|Circos parameters
@click.option('--window-size', type=int, default=1000000,
              show_default=True,
              help='Circos窗口大小|Circos window size (bp)')

@click.option('--aligner', default='minimap2',
              type=click.Choice(['minimap2', 'unimap']),
              show_default=True,
              help='比对工具|Aligner for homologous blocks')

# 高级参数|Advanced parameters
@click.option('--sg-assigned', default=None,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='已知亚基因组分配文件(跳过聚类)|Pre-assigned subgenome file (skip clustering)')

@click.option('--target', default=None,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='目标染色体文件|Target chromosomes file')

@click.option('--labels', multiple=True,
              help='基因组标签|Genome labels')

@click.option('--no-label', is_flag=True,
              help='不添加标签前缀|No label prefix')

@click.option('--custom-features', multiple=True,
              help='自定义特征FASTA|Custom feature FASTA files')

@click.option('--figfmt', default='pdf',
              type=click.Choice(['pdf', 'png']),
              show_default=True,
              help='图片格式|Figure format')

# 其他|Other
@click.option('--overwrite', is_flag=True,
              help='覆盖已有结果|Overwrite existing results')

@click.option('--cleanup', is_flag=True,
              help='清理临时文件|Clean up temporary files')

@click.option('--conda-env', default='SubPhaser',
              show_default=True,
              help='conda环境名称|Conda environment name')

def subphaser(genomes, nsg, sg_cfgs, parental_genomes, output_dir, prefix, threads,
              min_chrom_size, kmer_size, min_fold, min_freq, max_pval, replicates,
              test_method, disable_ltr, disable_circos, disable_blocks,
              just_core, ltr_detectors, mu, window_size, aligner,
              sg_assigned, target, labels, no_label, custom_features,
              figfmt, overwrite, cleanup, conda_env):
    """
    异源多倍体亚基因组分离|Phase subgenomes of allopolyploids

    基于重复k-mer的异源多倍体亚基因组自动分离和命名|
    Automatically phase and rename subgenomes based on repetitive kmers

    示例|Example: biopytools subphaser -i genome.fa --nsg 2 -o output/

    验证模式|Validation: biopytools subphaser -i genome.fa --nsg 2 --parental-genomes parentA.fa parentB.fa
    """

    subphaser_main = _lazy_import_subphaser_main()

    args = ['subphaser.py']

    # 必需参数|Required
    args.extend(['-i'] + list(genomes))
    args.extend(['--nsg', str(nsg)])

    # 输入模式|Input mode
    if sg_cfgs:
        args.extend(['-c'] + list(sg_cfgs))
    if parental_genomes:
        args.extend(['--parental-genomes'] + list(parental_genomes))

    # 输出|Output
    if output_dir != './subphaser_output':
        args.extend(['-o', output_dir])
    if prefix:
        args.extend(['--prefix', prefix])

    # 资源|Resources
    if threads != 24:
        args.extend(['-t', str(threads)])

    # 染色体过滤|Chromosome filtering
    if min_chrom_size != 1_000_000:
        args.extend(['--min-chrom-size', str(min_chrom_size)])

    # K-mer|K-mer
    if kmer_size != 15:
        args.extend(['-k', str(kmer_size)])
    if min_fold != 2.0:
        args.extend(['-f', str(min_fold)])
    if min_freq != 200:
        args.extend(['-q', str(min_freq)])

    # 聚类|Cluster
    if max_pval != 0.05:
        args.extend(['--max-pval', str(max_pval)])
    if replicates != 1000:
        args.extend(['--replicates', str(replicates)])
    if test_method != 'ttest_ind':
        args.extend(['--test-method', test_method])

    # 步骤控制|Step control
    if disable_ltr:
        args.append('--disable-ltr')
    if disable_circos:
        args.append('--disable-circos')
    if disable_blocks:
        args.append('--disable-blocks')
    if just_core:
        args.append('--just-core')

    # LTR
    if ltr_detectors:
        args.extend(['--ltr-detectors'] + list(ltr_detectors))
    if mu != 13e-9:
        args.extend(['--mu', str(mu)])

    # Circos
    if window_size != 1000000:
        args.extend(['--window-size', str(window_size)])
    if aligner != 'minimap2':
        args.extend(['--aligner', aligner])

    # 高级|Advanced
    if sg_assigned:
        args.extend(['--sg-assigned', sg_assigned])
    if target:
        args.extend(['--target', target])
    if labels:
        args.extend(['--labels'] + list(labels))
    if no_label:
        args.append('--no-label')
    if custom_features:
        args.extend(['--custom-features'] + list(custom_features))
    if figfmt != 'pdf':
        args.extend(['--figfmt', figfmt])

    # 其他|Other
    if overwrite:
        args.append('--overwrite')
    if cleanup:
        args.append('--cleanup')
    if conda_env != 'SubPhaser':
        args.extend(['--conda-env', conda_env])

    original_argv = sys.argv
    sys.argv = args

    try:
        subphaser_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
