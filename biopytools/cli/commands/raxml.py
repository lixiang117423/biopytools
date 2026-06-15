"""
RAxML系统发育分析命令|RAxML Phylogenetic Analysis Command
"""

import click
import sys
import os


def _lazy_import_raxml_main():
    """延迟加载raxml主函数|Lazy load raxml main function"""
    try:
        from ...raxml.main import main as raxml_main
        return raxml_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='RAxML系统发育分析工具|RAxML phylogenetic analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--sequence-file', '-s',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入序列文件(PHYLIP格式)|Input sequence file (PHYLIP format)')
@click.option('--output-name', '-n',
              required=True,
              type=str,
              help='输出文件名称|Output file name')
@click.option('--model', '-m',
              default='GTRGAMMA',
              type=str,
              show_default=True,
              help='替换模型|Substitution model (GTRGAMMA, PROTGAMMAWAG, etc.)')
@click.option('--categories', '-c',
              type=int,
              default=25,
              show_default=True,
              help='速率异质性类别数|Number of rate heterogeneity categories')
@click.option('--likelihood-epsilon', '-e',
              type=float,
              default=0.1,
              show_default=True,
              help='似然优化精度|Likelihood optimization precision')
@click.option('--algorithm', '-f',
              default='d',
              type=str,
              show_default=True,
              help='算法类型|Algorithm type (d=rapid hill-climbing, a=rapid bootstrap, etc.)')
@click.option('--parsimony-seed', '-p',
              type=int,
              help='简约法随机种子|Parsimony random seed')
@click.option('--runs', '-N',
              default='1',
              type=str,
              show_default=True,
              help='运行次数或bootstrap次数|Number of runs or bootstrap replicates')
@click.option('--bootstrap-seed', '-b',
              type=int,
              help='Bootstrap随机种子|Bootstrap random seed')
@click.option('--rapid-bootstrap-seed', '-x',
              type=int,
              help='快速bootstrap随机种子|Rapid bootstrap random seed')
@click.option('--bootstrap-convergence', '-I',
              type=click.Choice(['autoFC', 'autoMR', 'autoMRE', 'autoMRE_IGN']),
              help='Bootstrap收敛标准|Bootstrap convergence criterion')
@click.option('--bootstop-threshold', '-B',
              type=float,
              default=0.03,
              show_default=True,
              help='Bootstrap停止阈值|Bootstrap stop threshold')
@click.option('--bootstop-perms',
              type=int,
              default=100,
              show_default=True,
              help='Bootstrap停止检验次数|Bootstrap stop test permutations')
@click.option('--print-bootstrap-trees', '-k',
              is_flag=True,
              help='输出带分支长度的bootstrap树|Print bootstrap trees with branch lengths')
@click.option('--starting-tree', '-t',
              type=click.Path(exists=True),
              help='起始树文件|Starting tree file')
@click.option('--constraint-tree', '-g',
              type=click.Path(exists=True),
              help='约束树文件|Constraint tree file')
@click.option('--outgroup', '-o',
              type=str,
              help='外群名称(逗号分隔多个)|Outgroup name(s) (comma-separated)')
@click.option('--threads', '-T',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--memory-saving', '-U',
              is_flag=True,
              help='启用内存节省模式|Enable memory saving mode')
@click.option('--ml-search-convergence', '-D',
              is_flag=True,
              help='启用ML搜索收敛标准|Enable ML search convergence criterion')
@click.option('--random-starting-tree', '-d',
              is_flag=True,
              help='使用随机起始树|Use random starting tree')
@click.option('--disable-rate-heterogeneity', '-V',
              is_flag=True,
              help='禁用速率异质性模型|Disable rate heterogeneity model')
@click.option('--gamma-median', '-u',
              is_flag=True,
              help='使用GAMMA模型中位数|Use median for GAMMA model')
@click.option('--disable-pattern-compression', '-H',
              is_flag=True,
              help='禁用模式压缩|Disable pattern compression')
@click.option('--output-dir', '-w',
              default='./raxml_output',
              type=click.Path(),
              show_default=True,
              help='输出目录|Output directory')
@click.option('--raxml-path',
              default='raxmlHPC-PTHREADS',
              type=str,
              show_default=True,
              help='RAxML程序路径|RAxML program path')
@click.option('--no-seq-check',
              is_flag=True,
              help='跳过序列检查|Skip sequence checking')
@click.option('--silent',
              is_flag=True,
              help='静默模式|Silent mode')
def raxml(sequence_file, output_name, model, categories, likelihood_epsilon,
          algorithm, parsimony_seed, runs, bootstrap_seed, rapid_bootstrap_seed,
          bootstrap_convergence, bootstop_threshold, bootstop_perms, print_bootstrap_trees,
          starting_tree, constraint_tree, outgroup, threads, memory_saving,
          ml_search_convergence, random_starting_tree, disable_rate_heterogeneity,
          gamma_median, disable_pattern_compression, output_dir, raxml_path,
          no_seq_check, silent):
    """
    RAxML系统发育分析工具|RAxML Phylogenetic Analysis Tool

    基于RAxML进行系统发育树构建和bootstrap分析|Perform phylogenetic tree construction and bootstrap analysis using RAxML

    示例|Examples: biopytools raxml -s alignment.phy -n my_tree
    """

    # 延迟加载|Lazy load
    raxml_main = _lazy_import_raxml_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['raxml.py']

    # 必需参数|Required parameters
    args.extend(['-s', sequence_file])
    args.extend(['-n', output_name])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if model != 'GTRGAMMA':
        args.extend(['-m', model])

    if categories != 25:
        args.extend(['-c', str(categories)])

    if likelihood_epsilon != 0.1:
        args.extend(['-e', str(likelihood_epsilon)])

    if algorithm != 'd':
        args.extend(['-f', algorithm])

    if parsimony_seed is not None:
        args.extend(['-p', str(parsimony_seed)])

    if runs != '1':
        args.extend(['-N', runs])

    if bootstrap_seed is not None:
        args.extend(['-b', str(bootstrap_seed)])

    if rapid_bootstrap_seed is not None:
        args.extend(['-x', str(rapid_bootstrap_seed)])

    if bootstrap_convergence is not None:
        args.extend(['-I', bootstrap_convergence])

    if bootstop_threshold != 0.03:
        args.extend(['-B', str(bootstop_threshold)])

    if bootstop_perms != 100:
        args.extend(['--bootstop-perms', str(bootstop_perms)])

    if print_bootstrap_trees:
        args.append('-k')

    if starting_tree is not None:
        args.extend(['-t', starting_tree])

    if constraint_tree is not None:
        args.extend(['-g', constraint_tree])

    if outgroup is not None:
        args.extend(['-o', outgroup])

    if threads != 88:
        args.extend(['-T', str(threads)])

    if memory_saving:
        args.append('-U')

    if ml_search_convergence:
        args.append('-D')

    if random_starting_tree:
        args.append('-d')

    if disable_rate_heterogeneity:
        args.append('-V')

    if gamma_median:
        args.append('-u')

    if disable_pattern_compression:
        args.append('-H')

    if output_dir != './raxml_output':
        args.extend(['-w', output_dir])

    if raxml_path != 'raxmlHPC-PTHREADS':
        args.extend(['--raxml-path', raxml_path])

    if no_seq_check:
        args.append('--no-seq-check')

    if silent:
        args.append('--silent')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        raxml_main()
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
