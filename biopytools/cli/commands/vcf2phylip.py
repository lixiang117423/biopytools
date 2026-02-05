"""
VCF转换命令|VCF Converter Command
"""

import click
import sys


def _lazy_import_vcf2phylip_main():
    """延迟加载vcf2phylip主函数|Lazy load vcf2phylip main function"""
    try:
        from ...vcf2phylip.main import main as vcf2phylip_main
        return vcf2phylip_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


@click.command(
    short_help='VCF转换工具|VCF converter',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='输入VCF文件路径|Input VCF file path')
@click.option('--output', '-o',
              default='./converted_output',
              type=click.Path(),
              show_default=True,
              help='输出目录|Output directory')
@click.option('--output-prefix',
              type=str,
              help='输出文件名前缀|Output filename prefix')
@click.option('--min-samples-locus', '-m',
              default=4,
              type=int,
              show_default=True,
              help='位点最少样本数|Minimum samples per locus')
@click.option('--outgroup', '-g',
              default="",
              type=str,
              show_default=True,
              help='外群样本名称|Outgroup sample name')
@click.option('--phylip-disable', '-p',
              is_flag=True,
              help='禁用PHYLIP输出|Disable PHYLIP output')
@click.option('--fasta', '-f',
              is_flag=True,
              help='启用FASTA输出|Enable FASTA output')
@click.option('--nexus', '-n',
              is_flag=True,
              help='启用NEXUS输出|Enable NEXUS output')
@click.option('--nexus-binary', '-b',
              is_flag=True,
              help='启用二进制NEXUS输出|Enable binary NEXUS output')
@click.option('--resolve-IUPAC', '-r',
              is_flag=True,
              help='随机解析杂合子基因型|Resolve heterozygous genotypes')
@click.option('--write-used-sites', '-w',
              is_flag=True,
              help='保存筛选通过的位点坐标|Save used sites coordinates')
@click.option('--threads', '-t',
              default=12,
              type=int,
              show_default=True,
              help='线程数|Number of threads')
@click.version_option(version='v2.9.1', message='%(prog)s %(version)s')
def vcf2phylip(input, output, output_prefix, min_samples_locus, outgroup,
                  phylip_disable, fasta, nexus, nexus_binary, resolve_iupac,
                  write_used_sites, threads):
    """
    VCF转换工具|VCF Converter Tool

    将VCF格式SNP数据转换为PHYLIP/FASTA/NEXUS格式|Convert VCF SNPs to PHYLIP/FASTA/NEXUS formats

    示例|Examples:
    biopytools vcf2phylip -i variants.vcf -o converted_results
    """

    # 延迟加载|Lazy loading
    vcf2phylip_main = _lazy_import_vcf2phylip_main()

    # 构建参数列表|Build argument list
    args = ['vcf2phylip.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])

    # 输出参数|Output parameters
    if output != './converted_output':
        args.extend(['-o', output])

    if output_prefix:
        args.extend(['--output-prefix', output_prefix])

    # 转换参数|Conversion parameters
    if min_samples_locus != 4:
        args.extend(['-m', str(min_samples_locus)])

    if outgroup != "":
        args.extend(['-g', outgroup])

    # 输出格式控制|Output format control
    if phylip_disable:
        args.append('-p')

    if fasta:
        args.append('-f')

    if nexus:
        args.append('-n')

    if nexus_binary:
        args.append('-b')

    # 处理选项|Processing options
    if resolve_iupac:
        args.append('-r')

    if write_used_sites:
        args.append('-w')

    if threads != 88:
        args.extend(['-t', str(threads)])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        vcf2phylip_main()
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
