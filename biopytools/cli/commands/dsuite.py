"""
Dsuite命令|Dsuite Command
"""

import click
import sys
import os


def _lazy_import_dsuite_main():
    """延迟加载dsuite主函数|Lazy load dsuite main function"""
    try:
        from ...dsuite.main import main as dsuite_main
        return dsuite_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在|Validate file existence"""
    if _is_help_request():
        return file_path
    if not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(short_help='Dsuite D统计|Dsuite D-statistic analysis',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入VCF文件|Input VCF file path')
@click.option('--sets', '-s',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='SETS分组文件|SETS file path')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--prefix', '-p',
              default='dsuite',
              show_default=True,
              help='输出文件前缀|Output file prefix')
@click.option('--dsuite-bin',
              default='~/software/Dsuite/Build/Dsuite',
              show_default=True,
              help='Dsuite可执行文件路径|Dsuite binary path')
@click.option('--min-alleles',
              type=int,
              default=2,
              show_default=True,
              help='最小等位基因数|Min number of alleles')
@click.option('--max-alleles',
              type=int,
              default=2,
              show_default=True,
              help='最大等位基因数|Max number of alleles')
@click.option('--variant-type',
              default='snps',
              show_default=True,
              type=click.Choice(['snps', 'indels', 'both', 'none']),
              help='变异类型|Variant type')
@click.option('--bcftools',
              default='bcftools',
              show_default=True,
              help='bcftools命令路径|bcftools command path')
@click.option('--collect-stats',
              is_flag=True,
              default=False,
              help='收集VCF统计信息|Whether to collect VCF statistics')
def dsuite(input, sets, output_dir, prefix, dsuite_bin, min_alleles, max_alleles, variant_type, bcftools, collect_stats):
    """
    Dsuite D统计分析工具|Dsuite D-statistic Analysis Tool

    Dsuite D统计(ABBA-BABA测试)|Dsuite D-statistic (ABBA-BABA test)

    示例|Examples: biopytools dsuite -i variants.vcf.gz -s sets.txt -o output_dir
    """
    # 延迟加载|Lazy loading
    dsuite_main = _lazy_import_dsuite_main()

    # 构建参数列表|Build argument list
    args = ['dsuite.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-s', sets])
    args.extend(['-o', output_dir])

    # 可选参数|Optional parameters
    if prefix != 'dsuite':
        args.extend(['-p', prefix])

    if dsuite_bin != '~/software/Dsuite/Build/Dsuite':
        args.extend(['--dsuite-bin', dsuite_bin])

    if min_alleles != 2:
        args.extend(['--min-alleles', str(min_alleles)])

    if max_alleles != 2:
        args.extend(['--max-alleles', str(max_alleles)])

    if variant_type != 'snps':
        args.extend(['--variant-type', variant_type])

    if bcftools != 'bcftools':
        args.extend(['--bcftools', bcftools])

    if collect_stats:
        args.append('--collect-stats')

    # 保存和恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始main函数|Call original main function
        dsuite_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
