"""
VCF邻接树分析命令|VCF Neighbor-Joining Tree Analysis Command
"""

import click
import sys


def _lazy_import_vcf2nj_main():
    """延迟加载vcf2nj主函数|Lazy load vcf2nj main function"""
    try:
        from ...vcf2nj.main import main as vcf2nj_main
        return vcf2nj_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


@click.command(
    short_help='VCF邻接树分析工具|VCF neighbor-joining tree analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              type=click.Path(exists=True),
              help='输入VCF文件路径|Input VCF file path')
@click.option('--distance-matrix', '-d',
              type=click.Path(exists=True),
              help='已有距离矩阵文件路径|Existing distance matrix file path')
@click.option('--output', '-o',
              default='.',
              type=str,
              show_default=True,
              help='输出目录|Output directory')
@click.option('--prefix', '-p',
              default='vcf2nj',
              type=str,
              show_default=True,
              help='输出文件前缀|Output file prefix')
@click.option('--tree-output', '-t',
              type=click.Path(),
              help='系统发育树输出文件路径|Phylogenetic tree output file path')
@click.option('--outgroup',
              type=str,
              help='外群样本标签，多个用逗号分隔|Outgroup sample labels, comma-separated')
@click.option('--nw-reroot-path',
              default='~/miniforge3/envs/newick_utils_v.1.6/bin/nw_reroot',
              type=str,
              show_default=True,
              help='nw_reroot程序路径|nw_reroot program path')
@click.option('--vcf2dis-path',
              default='VCF2Dis',
              type=str,
              show_default=True,
              help='VCF2Dis程序路径|VCF2Dis program path')
@click.option('--working-dir', '-w',
              default='.',
              type=click.Path(),
              show_default=True,
              help='工作目录|Working directory')
@click.option('--skip-vcf2dis',
              is_flag=True,
              help='跳过VCF2Dis步骤|Skip VCF2Dis step')
def vcf2nj(input, distance_matrix, output, prefix, tree_output, outgroup, nw_reroot_path, vcf2dis_path, working_dir, skip_vcf2dis):
    """
    VCF邻接树分析工具|VCF Neighbor-Joining Tree Analysis Tool

    使用邻接算法从VCF文件进行系统发育分析|Perform phylogenetic analysis from VCF files using Neighbor-Joining algorithm

    示例|Examples: biopytools vcf2nj -i wild.snp.vcf -o output_dir -p wild_snp
    """

    # 参数验证|Parameter validation
    if skip_vcf2dis:
        if not distance_matrix:
            raise click.ClickException(
                "使用--skip-vcf2dis时必须指定--distance-matrix|"
                "--distance-matrix must be specified when using --skip-vcf2dis"
            )
    else:
        if not input:
            raise click.ClickException(
                "必须指定输入VCF文件|Input VCF file must be specified"
            )

    # 延迟加载|Lazy load
    vcf2nj_main = _lazy_import_vcf2nj_main()

    # 构建参数列表|Build argument list
    args = ['vcf2nj.py']

    # 输入参数|Input parameters
    if input:
        args.extend(['-i', input])

    if distance_matrix:
        args.extend(['-d', distance_matrix])

    # 输出参数|Output parameters
    if output != '.':
        args.extend(['-o', output])

    if prefix != 'vcf2nj':
        args.extend(['-p', prefix])

    if tree_output:
        args.extend(['-t', tree_output])

    # 重根化参数|Rerooting parameters
    if outgroup:
        args.extend(['--outgroup', outgroup])

    if nw_reroot_path != '~/miniforge3/envs/newick_utils_v.1.6/bin/nw_reroot':
        args.extend(['--nw-reroot-path', nw_reroot_path])

    # 工具路径|Tool paths
    if vcf2dis_path != 'VCF2Dis':
        args.extend(['--vcf2dis-path', vcf2dis_path])

    if working_dir != '.':
        args.extend(['-w', working_dir])

    # 执行选项|Execution options
    if skip_vcf2dis:
        args.append('--skip-vcf2dis')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主程序|Call main program
        vcf2nj_main()
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
