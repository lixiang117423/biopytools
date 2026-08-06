"""vcf2pav CLI 包装器|vcf2pav Click wrapper"""
import os
import sys

import click


def _lazy_import_main():
    """延迟导入主函数|Lazy import main"""
    try:
        from ...vcf2pav.main import main as vcf2pav_main
        return vcf2pav_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """是否帮助请求|Is help request"""
    return any(a in {"-h", "--help"} for a in sys.argv)


def _validate_path_exists(path):
    """校验路径存在|Validate path exists"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(short_help="VCF转PAV(Presence/Absence)矩阵|VCF to PAV matrix")
@click.option("-i", "--input", "input_vcf", required=True,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="输入 VCF 文件|Input VCF file")
@click.option("-o", "--output-dir", default="./vcf2pav_output",
              help="输出目录(默认./vcf2pav_output)|Output directory")
@click.option("--log-level", default="INFO",
              help="日志级别(默认INFO)|Log level (default INFO)")
def vcf2pav(**kwargs):
    """VCF转PAV(Presence/Absence Variation)矩阵,行=SV,列=样本,值=0/1.
    |Convert SURVIVOR-merged VCF to PAV matrix: rows=SVs, columns=samples, values=0/1.

    示例|Examples: biopytools vcf2pav -i pan_sv.survivor.vcf -o output_dir/
    """
    vcf2pav_main = _lazy_import_main()
    args = ["vcf2pav.py",
            "-i", kwargs["input_vcf"],
            "-o", kwargs["output_dir"],
            "--log-level", kwargs["log_level"]]
    original = sys.argv
    sys.argv = args
    try:
        rc = vcf2pav_main()
        sys.exit(rc if rc is not None else 0)
    except SystemExit as e:
        sys.exit(e.code if e.code is not None else 0)
    finally:
        sys.argv = original
