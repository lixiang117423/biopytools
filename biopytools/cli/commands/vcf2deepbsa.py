"""vcf2deepbsa CLI 包装器|vcf2deepbsa Click wrapper"""
import os
import sys

import click


def _lazy_import_main():
    """延迟导入主函数|Lazy import main"""
    try:
        from ...vcf2deepbsa.main import main as vcf2deepbsa_main
        return vcf2deepbsa_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """是否帮助请求|Is help request"""
    return any(a in {"-h", "--help"} for a in sys.argv)


def _validate_path_exists(path):
    """校验路径存在|Validate path exists"""
    if not _is_help_request() and path and not os.path.exists(os.path.expanduser(path)):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(short_help="VCF转DeepBSA输入CSV(提取AD)|VCF to DeepBSA input CSV (extract AD)")
@click.option("-i", "--input", "input_vcf", required=True,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="输入 VCF 文件(支持 .gz)|Input VCF file (.gz supported)")
@click.option("-o", "--output-dir", default="./vcf2deepbsa_output",
              help="输出目录(默认./vcf2deepbsa_output)|Output directory")
@click.option("--log-level", default="INFO",
              help="日志级别(默认INFO)|Log level (default INFO)")
def vcf2deepbsa(**kwargs):
    """VCF转DeepBSA输入CSV: 提取FORMAT中的AD(等位深度),输出无header的
    CHROM,POS,REF,ALT,AD对... 矩阵,供 deepbsa run/batch 直接使用.
    |Convert VCF to DeepBSA input CSV: extract AD (allele depth) from FORMAT
    into a headerless CHROM,POS,REF,ALT,AD-pairs... matrix.

    示例|Examples: biopytools vcf2deepbsa -i bsa_pools.vcf -o output_dir/
    """
    vcf2deepbsa_main = _lazy_import_main()
    args = ["vcf2deepbsa.py",
            "-i", kwargs["input_vcf"],
            "-o", kwargs["output_dir"],
            "--log-level", kwargs["log_level"]]
    original = sys.argv
    sys.argv = args
    try:
        rc = vcf2deepbsa_main()
        sys.exit(rc if rc is not None else 0)
    except SystemExit as e:
        sys.exit(e.code if e.code is not None else 0)
    finally:
        sys.argv = original
