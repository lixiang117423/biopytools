"""check_reads CLI 包装器|check_reads Click wrapper"""
import os
import sys

import click


def _lazy_import_main():
    """延迟导入主函数|Lazy import main"""
    try:
        from ...check_reads.main import main as check_reads_main
        return check_reads_main
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


@click.command(short_help="fastq完整性检查(gz/0字节/配对)|FASTQ integrity check "
                          "(gzip/empty/pairing)")
@click.option("-i", "--input", "input_dir", required=True,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="fastq 目录(逗号分隔多个,递归扫描)|FASTQ dir(s), "
                   "comma-separated, recursive")
@click.option("-o", "--output-dir", default="./check_reads_output",
              help="输出目录(默认./check_reads_output)|Output directory")
@click.option("-t", "--threads", default=12, type=int,
              help="并行线程数(默认12)|Parallel threads (default 12)")
@click.option("--log-level", default="INFO",
              help="日志级别(默认INFO)|Log level (default INFO)")
def check_reads(**kwargs):
    """fastq 完整性检查:gz 压缩完整性、0 字节文件、R1-R2 配对完整性
    |FASTQ integrity check: gzip integrity, empty files, R1-R2 pairing.

    示例|Examples: biopytools check-reads -i 2nd/clean/ -o check_out/ -t 88
    """
    check_reads_main = _lazy_import_main()
    args = ["check_reads.py",
            "-i", kwargs["input_dir"],
            "-o", kwargs["output_dir"],
            "-t", str(kwargs["threads"]),
            "--log-level", kwargs["log_level"]]
    original = sys.argv
    sys.argv = args
    try:
        rc = check_reads_main()
        sys.exit(rc if rc is not None else 0)
    except SystemExit as e:
        sys.exit(e.code if e.code is not None else 0)
    finally:
        sys.argv = original
