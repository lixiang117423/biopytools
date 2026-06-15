"""
IQ-TREE系统发育树分析CLI包装器|IQ-TREE Phylogenetic Tree Analysis CLI Wrapper
"""

import click
import sys
import os


def _lazy_import_iqtree_main():
    """延迟加载IQ-TREE主函数|Lazy load IQ-TREE main function"""
    try:
        from ...iqtree.main import main as iqtree_main
        return iqtree_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='IQ-TREE系统发育树分析|IQ-TREE Phylogenetic Tree Analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入比对文件|Input alignment file')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--prefix', '-p',
              required=True,
              help='输出文件前缀|Output file prefix')
@click.option('--model', '-m',
              help='进化模型|Evolutionary model')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--bootstrap', '-b',
              type=int,
              default=1000,
              show_default=True,
              help='Bootstrap重复次数|Bootstrap replicates')
@click.option('--boot-type',
              type=click.Choice(['ufboot', 'standard']),
              default='ufboot',
              show_default=True,
              help='Bootstrap类型|Bootstrap type')
@click.option('--save-boot-trees',
              is_flag=True,
              help='保存所有bootstrap树|Save all bootstrap trees')
@click.option('--outgroup',
              help='外群名称(逗号分隔)|Outgroup taxon names (comma-separated)')
@click.option('--constraint',
              type=click.Path(exists=True),
              help='约束树文件|Constraint tree file')
@click.option('--partition',
              type=click.Path(exists=True),
              help='分区文件|Partition file')
@click.option('--partition-mode',
              type=click.Choice(['edge-linked', 'edge-equal', 'edge-unlinked']),
              default='edge-linked',
              show_default=True,
              help='分区模式|Partition mode')
@click.option('--concordance',
              type=click.Path(exists=True),
              help='一致性因子基因树文件|Concordance factor gene tree file')
@click.option('--ancestral',
              is_flag=True,
              help='启用祖先状态重建|Enable ancestral state reconstruction')
@click.option('--seed',
              type=int,
              help='随机种子|Random seed')
@click.option('--runs',
              type=int,
              default=1,
              show_default=True,
              help='独立运行次数|Number of independent runs')
@click.option('--redo',
              is_flag=True,
              help='重新运行分析|Redo analysis')
@click.option('--iqtree-path',
              default='iqtree',
              show_default=True,
              help='IQ-TREE软件路径|IQ-TREE program path')
def iqtree(input, output, prefix, model, threads, bootstrap, boot_type, save_boot_trees,
           outgroup, constraint, partition, partition_mode,
           concordance, ancestral, seed, runs, redo, iqtree_path):
    """
    IQ-TREE系统发育树分析流程|IQ-TREE Phylogenetic Tree Analysis Pipeline

    基于多重序列比对构建系统发育树并进行分析
    Phylogenetic tree construction and analysis based on multiple sequence alignment

    示例|Examples: biopytools iqtree -i alignment.fasta -o tree_results -p my_tree
    """

    # 延迟加载|Lazy load
    iqtree_main = _lazy_import_iqtree_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['iqtree.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])
    args.extend(['-p', prefix])

    # 核心参数|Core parameters
    if model:
        args.extend(['-m', model])

    if threads != 88:
        args.extend(['-t', str(threads)])

    # Bootstrap参数|Bootstrap parameters
    if bootstrap != 1000:
        args.extend(['-b', str(bootstrap)])

    if boot_type != 'ufboot':
        args.extend(['--boot-type', boot_type])

    if save_boot_trees:
        args.append('--save-boot-trees')

    # 树选项|Tree options
    if outgroup:
        args.extend(['--outgroup', outgroup])

    if constraint:
        args.extend(['--constraint', constraint])

    # 分区分析|Partition analysis
    if partition:
        args.extend(['--partition', partition])

    if partition_mode != 'edge-linked':
        args.extend(['--partition-mode', partition_mode])

    # 高级功能|Advanced features
    if concordance:
        args.extend(['--concordance', concordance])

    if ancestral:
        args.append('--ancestral')

    # 其他参数|Other parameters
    if seed:
        args.extend(['--seed', str(seed)])

    if runs != 1:
        args.extend(['--runs', str(runs)])

    if redo:
        args.append('--redo')

    if iqtree_path != 'iqtree':
        args.extend(['--iqtree-path', iqtree_path])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        iqtree_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"意外错误|Unexpected error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
