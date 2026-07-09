"""
Phobius跨膜拓扑+信号肽预测命令|Phobius TM topology & signal peptide prediction command
"""

import click
import sys
import os


_DEFAULT_PHOBIUS_PATH = "~/miniforge3/envs/phobius_v.1.0.1/bin/phobius.pl"


def _lazy_import_phobius_main():
    """延迟加载phobius主函数|Lazy load phobius main function"""
    try:
        from ...phobius.main import main as phobius_main
        return phobius_main
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
    short_help="Phobius跨膜拓扑+信号肽预测|Phobius TM topology & signal peptide prediction",
    context_settings=dict(help_option_names=["-h", "--help"], max_content_width=120),
)
@click.option("--input", "-i",
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help="输入蛋白质FASTA|Input protein FASTA")
@click.option("--output-dir", "-o",
              required=True,
              type=click.Path(),
              help="输出目录|Output directory")
@click.option("--prefix",
              default=None,
              help="输出前缀(默认输入文件名)|Output prefix (default: input filename)")
@click.option("--phobius-path",
              default=_DEFAULT_PHOBIUS_PATH,
              show_default=True,
              help="phobius.pl路径|phobius.pl path")
def phobius(input, output_dir, prefix, phobius_path):
    """Phobius跨膜拓扑+信号肽预测|Phobius TM topology & signal peptide prediction

    示例|Examples: biopytools phobius -i proteins.fa -o output_dir/
    """
    phobius_main = _lazy_import_phobius_main()

    # 构建参数列表|Build argument list
    args = ["phobius.py"]
    args.extend(["-i", input])
    args.extend(["-o", output_dir])
    if prefix:
        args.extend(["--prefix", prefix])
    if phobius_path != _DEFAULT_PHOBIUS_PATH:
        args.extend(["--phobius-path", phobius_path])

    original_argv = sys.argv
    sys.argv = args
    try:
        phobius_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
