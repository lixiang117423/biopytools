"""vcf2splitstree CLI 包装器|vcf2splitstree Click wrapper"""
import os
import sys

import click


def _lazy_import_main():
    """延迟导入主函数|Lazy import main"""
    try:
        from ...vcf2splitstree.main import main as vcf2s_main
        return vcf2s_main
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


@click.command(short_help="VCF→SplitsTree6 距离矩阵|VCF to SplitsTree6 distance matrix")
@click.option("-i", "--input", "input_path", required=True,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="VCF 变异文件(.vcf/.vcf.gz)|VCF variants file")
@click.option("-o", "--output-dir", default="./vcf2splitstree_output",
              help="输出目录(默认./vcf2splitstree_output)|Output directory")
@click.option("--log-level", default="INFO",
              help="日志级别(默认INFO)|Log level (default INFO)")
def vcf2splitstree(input_path, output_dir, log_level):
    """VCF 变异文件 → SplitsTree6 可直接打开的距离矩阵 CSV
    |Convert VCF variants into a distance-matrix CSV that SplitsTree6 opens directly.

    在 SplitsTree6 GUI 中用 File→Open 打开输出的 CSV,会自动识别距离矩阵并
    运行 NeighborNet 网络分析。
    |Open the output CSV in SplitsTree6 GUI (File→Open); it auto-detects the
    distance matrix and runs NeighborNet.

    示例|Examples: biopytools vcf2splitstree -i variants.vcf.gz -o out/
    """
    vcf2s_main = _lazy_import_main()
    args = ["vcf2splitstree.py",
            "-i", input_path,
            "-o", output_dir,
            "--log-level", log_level]
    original = sys.argv
    sys.argv = args
    try:
        rc = vcf2s_main()
        sys.exit(rc if rc is not None else 0)
    except SystemExit as e:
        sys.exit(e.code if e.code is not None else 0)
    finally:
        sys.argv = original
