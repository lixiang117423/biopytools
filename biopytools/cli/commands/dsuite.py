"""
Dsuite Command
Dsuite D统计量分析命令
"""

import click
import sys
import os


def _lazy_import_dsuite_main():
    """懒加载dsuite main函数 | Lazy load dsuite main function"""
    try:
        from ...dsuite.main import main as dsuite_main
        return dsuite_main
    except ImportError as e:
        click.echo(f"[ERROR] 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path


@click.command(short_help="Dsuite D统计量分析工具",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--vcf', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='[FILE] 输入VCF文件路径 | Input VCF file path')
@click.option('--sets', '-s',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='[FILE] SETS分组文件路径 | SETS file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='[DIR] 输出目录 | Output directory')
@click.option('--prefix', '-p',
              default='dsuite',
              help='[STR] 输出文件前缀 (默认: dsuite) | Output file prefix (default: dsuite)')
@click.option('--dsuite-bin',
              default='/share/org/YZWL/yzwl_lixg/software/Dsuite/Build/Dsuite',
              help='[FILE] Dsuite可执行文件路径 | Dsuite binary path')
@click.option('--min-alleles',
              type=int,
              default=2,
              help='[INT] 最小等位基因数 (默认: 2) | Min number of alleles (default: 2)')
@click.option('--max-alleles',
              type=int,
              default=2,
              help='[INT] 最大等位基因数 (默认: 2) | Max number of alleles (default: 2)')
@click.option('--variant-type',
              default='snps',
              type=click.Choice(['snps', 'indels', 'both', 'none']),
              help='[STR] 变异类型 (默认: snps) | Variant type (default: snps)')
@click.option('--bcftools',
              default='bcftools',
              help='[CMD] bcftools命令路径 | bcftools command path')
def dsuite(vcf, sets, output, prefix, dsuite_bin, min_alleles, max_alleles, variant_type, bcftools):
    """
    Dsuite D统计量分析工具

    使用Dsuite进行基因渗入分析，计算D统计量(ABBA-BABA test)
    检测群体间的基因渗入事件。

    示例 | Examples:

    基本用法:
      biopytools dsuite -i variants.vcf.gz -s sets.txt -o output_dir

    指定输出前缀:
      biopytools dsuite -i variants.vcf.gz -s sets.txt -o output_dir -p analysis

    自定义过滤参数:
      biopytools dsuite -i variants.vcf.gz -s sets.txt -o output_dir \\
          --min-alleles 2 --max-alleles 2 --variant-type snps

    输出说明 | Output Description:

    工具会生成以下文件:
      - dsuite_BBAA.txt: D统计量和f4-ratio结果
      - dsuite_Dmin.txt: D最小值结果
      - dsuite_tree.txt: 按树结构排列的结果
      - dsuite_analysis_*.log: 运行日志

    SETS文件格式 | SETS File Format:

    SETS.txt文件定义样本所属的群体，格式为:
      SAMPLE_ID    POPULATION_ID

    例如:
      Ind1    Species1
      Ind2    Species1
      Ind3    Species2
      Ind4    Species2
      Ind5    Outgroup

    注意事项 | Notes:

      - 至少需要一个样本标记为Outgroup
      - 默认只分析双等位SNP
      - 使用bcftools进行VCF过滤
      - 分析结果可用于检测基因渗入事件
    """
    # 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    dsuite_main = _lazy_import_dsuite_main()

    # 构建参数列表传递给原始main函数 | Build argument list for original main function
    args = ['dsuite.py']

    # 必需参数 | Required parameters
    args.extend(['-i', vcf])
    args.extend(['-s', sets])
    args.extend(['-o', output])

    # 可选参数（只在非默认值时添加，减少命令行长度）| Optional parameters (add only when non-default)
    if prefix != 'dsuite':
        args.extend(['-p', prefix])

    if dsuite_bin != '/share/org/YZWL/yzwl_lixg/software/Dsuite/Build/Dsuite':
        args.extend(['--dsuite-bin', dsuite_bin])

    if min_alleles != 2:
        args.extend(['--min-alleles', str(min_alleles)])

    if max_alleles != 2:
        args.extend(['--max-alleles', str(max_alleles)])

    if variant_type != 'snps':
        args.extend(['--variant-type', variant_type])

    if bcftools != 'bcftools':
        args.extend(['--bcftools', bcftools])

    # 保存并恢复sys.argv | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 | Call original main function
        dsuite_main()
    except SystemExit as e:
        # 处理程序正常退出 | Handle normal program exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n[WARNING] 用户中断操作 | User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"[ERROR] 运行错误 | Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
