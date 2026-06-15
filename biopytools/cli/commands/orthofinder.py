"""
OrthoFinder泛基因组分析命令|OrthoFinder Pangenome Analysis Command
"""

import click
import sys
import os


def _lazy_import_orthofinder_main():
    """延迟加载orthofinder主函数|Lazy load orthofinder main function"""
    try:
        from ...orthofinder.main import main as orthofinder_main
        return orthofinder_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_directory_exists(dir_path):
    """验证目录存在(仅在非帮助模式)|Validate directory existence (only in non-help mode)"""
    if _is_help_request():
        return dir_path
    if not os.path.exists(dir_path):
        raise click.BadParameter(f"目录不存在|Directory does not exist: {dir_path}")
    if not os.path.isdir(dir_path):
        raise click.BadParameter(f"路径不是目录|Path is not a directory: {dir_path}")
    return dir_path


@click.command(
    short_help='OrthoFinder泛基因组分析工具|OrthoFinder pangenome analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              help='输入蛋白质序列文件目录|Input protein sequence files directory')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--project-name', '-n',
              type=str,
              help='项目名称|Project name')
@click.option('--softcore-threshold',
              type=int,
              default=1,
              show_default=True,
              help='Softcore基因缺失阈值|Softcore gene missing threshold')
@click.option('--dispensable-threshold',
              type=int,
              default=1,
              show_default=True,
              help='Dispensable基因缺失阈值|Dispensable gene missing threshold')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--search', '--search-program', '-s',
              type=click.Choice(['blast', 'diamond', 'diamond_ultra_sens', 'mmseqs', 'blast_nucl']),
              default='blast',
              show_default=True,
              help='序列搜索程序|Sequence search program')
@click.option('--mcl-inflation',
              type=float,
              default=1.2,
              show_default=True,
              help='MCL膨胀参数|MCL inflation parameter')
@click.option('--dna', '-d',
              is_flag=True,
              help='输入序列为DNA|Input sequences are DNA')
@click.option('--basic-only',
              is_flag=True,
              help='仅执行基础分析|Perform basic analysis only')
@click.option('--generate-trees',
              is_flag=True,
              help='生成系统发育树|Generate phylogenetic trees')
@click.option('--msa-program',
              type=click.Choice(['mafft', 'muscle']),
              default='mafft',
              show_default=True,
              help='多序列比对程序|Multiple sequence alignment program')
@click.option('--tree-program',
              type=click.Choice(['fasttree', 'fasttree_fastest', 'raxml', 'raxml-ng', 'iqtree']),
              default='fasttree',
              show_default=True,
              help='系统发育树推断程序|Phylogenetic tree inference program')
@click.option('--orthofinder-path',
              default='orthofinder',
              show_default=True,
              help='OrthoFinder程序路径|OrthoFinder program path')
@click.option('--force',
              is_flag=True,
              help='强制重新分析覆盖已有结果|Force reanalysis overwriting existing results')
@click.option('--skip-orthofinder',
              is_flag=True,
              help='跳过OrthoFinder步骤直接分类|Skip OrthoFinder step and go directly to classification')
@click.option('--disable-rarefaction',
              is_flag=True,
              help='禁用稀释曲线分析|Disable rarefaction curve analysis')
@click.option('--disable-single-copy',
              is_flag=True,
              help='禁用单拷贝基因分析|Disable single copy gene analysis')
@click.option('--no-plots',
              is_flag=True,
              help='不生成图表|Do not generate plots')
@click.option('--rarefaction-iterations',
              type=int,
              default=100,
              show_default=True,
              help='稀释分析迭代次数|Rarefaction analysis iterations')
@click.option('--single-copy-format',
              type=click.Choice(['by_orthogroup', 'by_genome', 'both']),
              default='both',
              show_default=True,
              help='单拷贝基因序列输出格式|Single copy gene sequence output format')
@click.option('--plot-format',
              type=click.Choice(['png', 'pdf', 'svg']),
              default='png',
              show_default=True,
              help='图表格式|Plot format')
def orthofinder(input, output, project_name, softcore_threshold, dispensable_threshold,
                threads, search, mcl_inflation, dna, basic_only, generate_trees,
                msa_program, tree_program, orthofinder_path, force, skip_orthofinder,
                disable_rarefaction, disable_single_copy, no_plots,
                rarefaction_iterations, single_copy_format, plot_format):
    """
    OrthoFinder泛基因组分析工具|OrthoFinder Pangenome Analysis Tool

    基于OrthoFinder进行泛基因组分析和直系同源群推断
    Perform pangenome analysis and orthogroup inference using OrthoFinder

    示例|Examples: biopytools orthofinder -i protein_sequences/ -o pangenome_results
    """

    # 延迟加载|Lazy load
    orthofinder_main = _lazy_import_orthofinder_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['orthofinder.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if project_name:
        args.extend(['-n', project_name])

    if softcore_threshold != 1:
        args.extend(['--softcore-threshold', str(softcore_threshold)])

    if dispensable_threshold != 1:
        args.extend(['--dispensable-threshold', str(dispensable_threshold)])

    if threads != 12:
        args.extend(['-t', str(threads)])

    if search != 'blast':
        args.extend(['-s', search])

    if mcl_inflation != 1.2:
        args.extend(['--mcl-inflation', str(mcl_inflation)])

    if msa_program != 'mafft':
        args.extend(['--msa-program', msa_program])

    if tree_program != 'fasttree':
        args.extend(['--tree-program', tree_program])

    if orthofinder_path != 'orthofinder':
        args.extend(['--orthofinder-path', orthofinder_path])

    if force:
        args.append('--force')

    if skip_orthofinder:
        args.append('--skip-orthofinder')

    if disable_rarefaction:
        args.append('--disable-rarefaction')

    if disable_single_copy:
        args.append('--disable-single-copy')

    if no_plots:
        args.append('--no-plots')

    if rarefaction_iterations != 100:
        args.extend(['--rarefaction-iterations', str(rarefaction_iterations)])

    if single_copy_format != 'both':
        args.extend(['--single-copy-format', single_copy_format])

    if plot_format != 'png':
        args.extend(['--plot-format', plot_format])

    # 布尔选项|Boolean options
    if dna:
        args.append('-d')

    if generate_trees:
        args.append('--generate-trees')

    if basic_only:
        args.append('--basic-only')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        orthofinder_main()
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
