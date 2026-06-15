"""
系统发育树样品选择命令|Phylogenetic Tree Sample Selection Command

基于层级关系和PCA的智能选择算法|Intelligent selection based on hierarchy and PCA
"""

import click
import sys


def _lazy_import_phylo_selector_main():
    """延迟加载phylo_selector主函数|Lazy load phylo_selector main function"""
    try:
        from ...phylo_selector.main import main as phylo_selector_main
        return phylo_selector_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if file_path:
        import os
        if not os.path.exists(file_path):
            raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='系统发育树样品选择（基于层级+PCA）|Phylogenetic tree sample selection (Hierarchy+PCA)',
    context_settings=dict(help_option_names=['--help'], max_content_width=120)
)
@click.option('-x', '--hierarchy-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='层级关系文件（必需）|Hierarchy file (required)')
@click.option('-p', '--pca-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='PCA坐标文件（必需）|PCA coordinates file (required)')
@click.option('-o', '--output-prefix',
              required=True,
              help='输出文件前缀|Output file prefix')
@click.option('-n', '--n-samples',
              type=int,
              default=150,
              show_default=True,
              help='选择样品总数|Total number of samples to select')
@click.option('-g', '--group-file',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='样品分组表文件|Sample group table file')
@click.option('--hierarchy-level',
              type=int,
              default=10,
              show_default=True,
              help='层级深度（用于分组）|Hierarchy level for grouping')
@click.option('--pca-threshold',
              type=float,
              default=0.0001,
              show_default=True,
              help='PCA去重阈值|PCA dedup threshold')
@click.option('--min-samples-per-group',
              type=int,
              default=1,
              show_default=True,
              help='每组最小样品数|Minimum samples per group')
@click.option('--no-report',
              is_flag=True,
              help='不生成详细报告|Do not generate detailed report')
@click.option('--no-csv',
              is_flag=True,
              help='不生成CSV文件|Do not generate CSV file')
@click.option('--no-visualization',
              is_flag=True,
              help='不生成可视化|Do not generate visualization')
def phylo_selector(hierarchy_file, pca_file, output_prefix, n_samples,
                   group_file, hierarchy_level, pca_threshold,
                   min_samples_per_group, no_report, no_csv, no_visualization):
    """
    系统发育树样品选择工具（基于层级关系和PCA）|Phylogenetic Tree Sample Selector (Hierarchy+PCA)

    从层级关系文件中选择代表性样品，使用PCA去重确保多样性
    Select representative samples from hierarchy file, using PCA dedup to ensure diversity

    示例|Examples: biopytools phylo-selector -x hierarchy.xlsx -p pca.txt -o selected_samples -n 150
    """

    # 延迟加载|Lazy loading
    phylo_selector_main = _lazy_import_phylo_selector_main()

    # 构建参数列表|Build argument list
    args = ['phylo_selector.py']

    # 必需参数|Required parameters
    args.extend(['--hierarchy-file', hierarchy_file])
    args.extend(['--pca-file', pca_file])
    args.extend(['--output-prefix', output_prefix])

    # 可选参数|Optional parameters
    if n_samples != 150:
        args.extend(['--n-samples', str(n_samples)])

    if group_file:
        args.extend(['--group-file', group_file])

    if hierarchy_level != 10:
        args.extend(['--hierarchy-level', str(hierarchy_level)])

    if pca_threshold != 0.0001:
        args.extend(['--pca-threshold', str(pca_threshold)])

    if min_samples_per_group != 1:
        args.extend(['--min-samples-per-group', str(min_samples_per_group)])

    # 布尔选项|Boolean options
    if no_report:
        args.append('--no-report')

    if no_csv:
        args.append('--no-csv')

    if no_visualization:
        args.append('--no-visualization')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        phylo_selector_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
