"""
VCF2Gene变异基因注释|VCF2Gene Variant Gene Annotation Command

"""

import click
import sys
import os


def _lazy_import_vcf2gene_main():
    """延迟加载vcf2gene主函数|Lazy load vcf2gene main function"""
    try:
        from ...vcf2gene.main import main as vcf2gene_main
        return vcf2gene_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='VCF变异基因注释工具|VCF variant gene annotation tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--vcf', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入VCF文件|Input VCF file path')
@click.option('--gff', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入GFF注释文件|Input GFF annotation file path')
@click.option('--output', '-o',
              required=True,
              help='输出结果文件|Output result file path')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads (reserved for future use)')
def vcf2gene(vcf, gff, output, threads):
    """
    VCF变异基因注释工具|VCF Variant Gene Annotation Tool

    基于GFF注释文件，标注VCF变异位点所在的基因区域（外显子/内含子/基因间区）|
    Annotate VCF variant positions with gene regions (exon/intron/intergenic) based on GFF annotation

    示例|Example: biopytools vcf2gene -i variants.vcf -g annotation.gff -o annotated.txt
    """

    # 延迟加载|Lazy loading
    vcf2gene_main = _lazy_import_vcf2gene_main()

    # 构建参数列表|Build argument list
    args = ['vcf2gene.py']

    # 必需参数|Required parameters
    args.extend(['-i', vcf])
    args.extend(['-g', gff])
    args.extend(['-o', output])

    # 可选参数|Optional parameters
    if threads != 12:
        args.extend(['-t', str(threads)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        vcf2gene_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
