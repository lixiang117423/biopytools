"""
RAxML-NG系统发育树分析CLI包装器|RAxML-NG Phylogenetic Tree Analysis CLI Wrapper
"""

import click
import sys
import os


def _lazy_import_raxml_ng_main():
    """延迟加载RAxML-NG主函数|Lazy load RAxML-NG main function"""
    try:
        from ...raxml_ng.main import main as raxml_ng_main
        return raxml_ng_main
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
    short_help='RAxML-NG系统发育树分析|RAxML-NG Phylogenetic Tree Analysis',
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
              default=None,
              help='输出文件前缀(默认输入文件名)|Output prefix (default: input filename)')
@click.option('--mode',
              type=click.Choice(['all', 'search', 'support']),
              default='all',
              show_default=True,
              help='分析模式|Analysis mode')
@click.option('--model', '-m',
              help='进化模型(不指定则自适应)|Evolutionary model (auto if not specified)')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--bs-trees', '-b',
              default='1000',
              show_default=True,
              help='Bootstrap重复数或autoMRE{N}|Bootstrap replicates or autoMRE{N}')
@click.option('--bs-metric',
              default='fbp',
              show_default=True,
              help='支持值类型: fbp/tbe/sh/ebg/rbs|Branch support metric')
@click.option('--tree',
              type=click.Path(exists=True),
              help='起始树(search)/参考树(support)|Starting tree (search) / reference tree (support)')
@click.option('--outgroup',
              help='外群名称(逗号分隔)|Outgroup taxon names (comma-separated)')
@click.option('--seed',
              type=int,
              help='随机种子|Random seed')
@click.option('--redo',
              is_flag=True,
              help='覆盖已有结果,忽略checkpoint|Overwrite results, ignore checkpoint')
@click.option('--raxml-ng-path',
              help='RAxML-NG软件路径|RAxML-NG program path')
def raxml_ng(input, output, prefix, mode, model, threads, bs_trees, bs_metric,
             tree, outgroup, seed, redo, raxml_ng_path):
    """
    RAxML-NG系统发育树分析流程|RAxML-NG Phylogenetic Tree Analysis Pipeline

    基于多重序列比对用最大似然法构建系统发育树
    Maximum likelihood phylogenetic tree construction from multiple sequence alignment

    示例|Examples: biopytools raxml-ng -i alignment.fasta -o tree_results -p my_tree
    """

    raxml_ng_main = _lazy_import_raxml_ng_main()

    args = ['raxml_ng.py']

    # 必需参数|Required
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 核心参数(仅非默认值透传,保持 help 简洁)|Core (only non-defaults)
    if prefix:
        args.extend(['-p', prefix])
    if mode != 'all':
        args.extend(['--mode', mode])
    if model:
        args.extend(['-m', model])
    if threads != 12:
        args.extend(['-t', str(threads)])

    # Bootstrap|Bootstrap
    if bs_trees != '1000':
        args.extend(['-b', bs_trees])
    if bs_metric != 'fbp':
        args.extend(['--bs-metric', bs_metric])

    # 树选项|Tree options
    if tree:
        args.extend(['--tree', tree])
    if outgroup:
        args.extend(['--outgroup', outgroup])

    # 其他|Other
    if seed is not None:
        args.extend(['--seed', str(seed)])
    if redo:
        args.append('--redo')
    if raxml_ng_path:
        args.extend(['--raxml-ng-path', raxml_ng_path])

    original_argv = sys.argv
    sys.argv = args

    try:
        raxml_ng_main()
    except SystemExit as e:
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
