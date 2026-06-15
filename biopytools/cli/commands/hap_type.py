"""
单倍型提取命令|Haplotype Extraction Command
"""

import click
import sys
from pathlib import Path


def _lazy_import_hap_type_main():
    """延迟加载hap_type主函数|Lazy load hap_type main function"""
    try:
        from ...hap_type.main import main as hap_type_main
        return hap_type_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not Path(file_path).exists():
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='单倍型提取工具|Haplotype extraction tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--vcf',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='VCF变异文件|VCF variant file')
@click.option('-r', '--region',
              required=True,
              help='基因组区间(chr:start-end)或BED文件|Genomic interval or BED file')
@click.option('-o', '--output',
              default=None,
              help='输出文件前缀(可选，默认自动生成)|Output file prefix (auto-generated if omitted)')
@click.option('--hetero-remove/--no-hetero-remove',
              default=False,
              show_default=True,
              help='去除杂合位点|Remove heterozygous sites')
@click.option('--na-drop/--no-na-drop',
              default=True,
              show_default=True,
              help='去除缺失位点|Remove missing sites')
@click.option('--hap-prefix',
              default='H',
              show_default=True,
              help='单倍型ID前缀|Haplotype ID prefix')
@click.option('--pad',
              type=int,
              default=3,
              show_default=True,
              help='单倍型ID位数|Haplotype ID padding digits')
def hap_type(vcf, region, output, hetero_remove, na_drop, hap_prefix, pad):
    """
    单倍型提取工具|Haplotype Extraction Tool

    从VCF文件指定区间提取单倍型，输出兼容geneHapR格式|Extract haplotypes from VCF in specified region, geneHapR-compatible output

    通过 -r 指定单个区间(chr:start-end)或BED文件批量处理|Use -r for single region or BED file for batch processing

    示例|Example:
    biopytools hap-type -i sample.vcf -r regions.bed -o result
    """

    # BED文件存在性校验|Validate BED file existence
    import re
    if region and not re.match(r'^[^:]+:\s*\d+\s*[-]\s*\d+$', region.strip()):
        _validate_file_exists(region)

    hap_type_main = _lazy_import_hap_type_main()

    class Args:
        def __init__(self):
            self.vcf = vcf
            self.region = region
            self.output = output or ""
            self.hetero_remove = hetero_remove
            self.na_drop = na_drop
            self.hap_prefix = hap_prefix
            self.pad = pad

    args = Args()

    try:
        exit_code = hap_type_main(args)
        sys.exit(exit_code)
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)
