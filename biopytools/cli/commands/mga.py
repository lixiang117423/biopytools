"""MGA共识基因组组装命令|MGA consensus genome assembly command"""

import sys
import os

import click


_DEFAULT_MGA_PATH = "~/software/MGA/consensusLJA/bin/MGA"


def _lazy_import_mga_main():
    """延迟加载mga主函数|Lazy load mga main function"""
    try:
        from ...mga.main import main as mga_main
        return mga_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否帮助请求|Check if help request"""
    help_flags = {"-h", "--help"}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(非帮助模式)|Validate file exists (non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help="MGA共识基因组组装(HiFi)|MGA consensus genome assembly (HiFi)",
    context_settings=dict(help_option_names=["-h", "--help"], max_content_width=120),
)
@click.option("-r", "--reads",
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help="HiFi reads(fasta/fastq,可gz)|HiFi reads (fasta/fastq, may be gz)")
@click.option("-o", "--output-dir",
              required=True,
              type=click.Path(),
              help="输出目录|Output directory")
@click.option("-t", "--threads",
              default=50, show_default=True,
              help="线程数|Threads")
@click.option("--mga-path",
              default=_DEFAULT_MGA_PATH, show_default=True,
              help="MGA二进制路径|MGA binary path")
@click.option("--conda-env",
              default="mga", show_default=True,
              help="conda环境名|conda env name")
@click.option("--dry-run", is_flag=True, default=False,
              help="只打印命令不执行|Print command without executing")
def mga(reads, output_dir, threads, mga_path, conda_env, dry_run):
    """MGA共识基因组组装|MGA consensus genome assembly

    示例|Examples: biopytools mga -r reads.fastq.gz -o out_dir/
    """
    mga_main = _lazy_import_mga_main()

    # 构建参数列表|Build argument list
    args = ["mga.py"]
    args.extend(["-r", reads])
    args.extend(["-o", output_dir])
    args.extend(["-t", str(threads)])
    if mga_path != _DEFAULT_MGA_PATH:
        args.extend(["--mga-path", mga_path])
    if conda_env != "mga":
        args.extend(["--conda-env", conda_env])
    if dry_run:
        args.append("--dry-run")

    original_argv = sys.argv
    sys.argv = args
    try:
        mga_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
