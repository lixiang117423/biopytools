"""
VCF转树命令|VCF to Tree Command
"""

import click
import sys
import os


def _lazy_import_vcf2tree_main():
    """延迟加载vcf2tree主函数|Lazy load vcf2tree main function"""
    try:
        from ...vcf2tree.main import main as vcf2tree_main
        return vcf2tree_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在|Validate file existence"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(
            f"文件不存在|File does not exist: {file_path}"
        )
    return file_path


@click.command(
    short_help='VCF转系统发育树|VCF to phylogenetic tree',
    context_settings=dict(
        help_option_names=['-h', '--help'], max_content_width=120
    )
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: (
                  _validate_file_exists(value) if value else None
              ),
              help='输入VCF文件路径|Input VCF file path')
@click.option('--method', '-m',
              type=click.Choice(['fasttree', 'iqtree']),
              default='iqtree',
              show_default=True,
              help='建树方法|Tree method (默认: iqtree|default: iqtree)')
@click.option('--output-dir', '-o',
              default='./vcf2tree_output',
              type=click.Path(),
              show_default=True,
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--outgroup', '-g',
              default='',
              show_default=True,
              help='外群样本名称|Outgroup sample name')
@click.option('--min-samples-locus',
              type=int,
              default=4,
              show_default=True,
              help='位点最少样本数|Minimum samples per locus')
@click.option('--fasttree-path',
              default='~/.local/bin/FastTree',
              show_default=True,
              help='FastTree软件路径|FastTree software path')
@click.option('--fasttree-params',
              default='',
              show_default=True,
              help='FastTree额外参数|Additional FastTree parameters')
@click.option('--iqtree-path',
              default=None,
              help='IQ-TREE软件路径|IQ-TREE software path')
@click.option('--iqtree-bootstrap',
              type=int,
              default=1000,
              show_default=True,
              help='IQ-TREE UFBoot重复次数|IQ-TREE UFBoot replicates')
@click.option('--iqtree-model',
              default=None,
              help='IQ-TREE模型(默认ModelFinder)|IQ-TREE model (default: ModelFinder)')
@click.option('--no-asc',
              is_flag=True,
              default=False,
              help='关闭SNP数据的ASC校正(默认开启)|Disable ASC correction for SNP data (on by default)')
def vcf2tree(input, method, output_dir, threads, outgroup, min_samples_locus,
             fasttree_path, fasttree_params, iqtree_path,
             iqtree_bootstrap, iqtree_model, no_asc):
    """
    VCF转系统发育树工具|VCF to Phylogenetic Tree Tool

    从VCF SNP数据直接构建系统发育树，默认使用IQ-TREE (ModelFinder + UFBoot)
    |Direct phylogenetic tree from VCF SNPs, default IQ-TREE (ModelFinder + UFBoot)

    示例|Examples: biopytools vcf2tree -i variants.vcf -o tree_output
    """

    vcf2tree_main = _lazy_import_vcf2tree_main()

    args = ['vcf2tree.py']
    args.extend(['-i', input])

    if method != 'iqtree':
        args.extend(['--method', method])
    if output_dir != './vcf2tree_output':
        args.extend(['-o', output_dir])
    if threads != 12:
        args.extend(['-t', str(threads)])
    if outgroup:
        args.extend(['-g', outgroup])
    if min_samples_locus != 4:
        args.extend(['--min-samples-locus', str(min_samples_locus)])
    if fasttree_path != '~/.local/bin/FastTree':
        args.extend(['--fasttree-path', fasttree_path])
    if fasttree_params:
        args.extend(['--fasttree-params', fasttree_params])
    if iqtree_path:
        args.extend(['--iqtree-path', iqtree_path])
    if iqtree_bootstrap != 1000:
        args.extend(['--iqtree-bootstrap', str(iqtree_bootstrap)])
    if iqtree_model:
        args.extend(['--iqtree-model', iqtree_model])
    if no_asc:
        args.extend(['--no-asc'])

    original_argv = sys.argv
    sys.argv = args

    try:
        vcf2tree_main()
    except SystemExit as e:
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
