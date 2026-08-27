"""phy2fa CLI 包装器|phy2fa Click wrapper"""
import os
import sys

import click


def _lazy_import_main():
    """延迟导入主函数|Lazy import main"""
    try:
        from ...phy2fa.main import main as phy2fa_main
        return phy2fa_main
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


@click.command(short_help="Phylip→FASTA 转换(sequential/interleaved)|"
                          "Phylip to FASTA conversion (sequential/interleaved)")
@click.option("-i", "--input", "input_path", required=True,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="Phylip 文件(.phy/.phylip/.gz)|Phylip file")
@click.option("-o", "--output-dir", default="./phy2fa_output",
              help="输出目录(默认./phy2fa_output)|Output directory")
@click.option("--line-width", default=60, type=int,
              help="FASTA 换行宽度,0=不换行(默认60)|FASTA line wrap width, "
                   "0=no wrap (default 60)")
@click.option("--log-level", default="INFO",
              help="日志级别(默认INFO)|Log level (default INFO)")
def phy2fa(input_path, output_dir, line_width, log_level):
    """Phylip 序列矩阵 → FASTA,自动识别 sequential/interleaved 两种排版
    |Convert a Phylip sequence matrix to FASTA; auto-detects both layouts.

    示例|Examples: biopytools phy2fa -i aln.phy -o fa_out/
    """
    phy2fa_main = _lazy_import_main()
    args = ["phy2fa.py",
            "-i", input_path,
            "-o", output_dir,
            "--line-width", str(line_width),
            "--log-level", log_level]
    original = sys.argv
    sys.argv = args
    try:
        rc = phy2fa_main()
        sys.exit(rc if rc is not None else 0)
    except SystemExit as e:
        sys.exit(e.code if e.code is not None else 0)
    finally:
        sys.argv = original
