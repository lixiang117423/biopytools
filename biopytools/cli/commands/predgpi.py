"""
PredGPI GPI锚定蛋白预测命令|PredGPI GPI-anchored protein prediction command
"""

import click
import sys
import os


_DEFAULT_PREDGPI_HOME = "~/software/predgpi"


def _lazy_import_predgpi_main():
    """延迟加载predgpi主函数|Lazy load predgpi main function"""
    try:
        from ...predgpi.main import main as predgpi_main
        return predgpi_main
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
    short_help="PredGPI GPI锚定蛋白预测|PredGPI GPI-anchor prediction",
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
@click.option("--predgpi-home",
              default=_DEFAULT_PREDGPI_HOME,
              show_default=True,
              help="predgpi安装目录|predgpi install directory")
@click.option("--conservative",
              is_flag=True,
              help="使用保守omega模型|Use conservative omega model")
@click.option("--prefix",
              default=None,
              help="输出前缀(默认输入文件名)|Output prefix (default: input filename)")
def predgpi(input, output_dir, predgpi_home, conservative, prefix):
    """PredGPI GPI锚定蛋白预测|PredGPI GPI-anchor prediction

    示例|Examples: biopytools predgpi -i proteins.fa -o output_dir/
    """
    predgpi_main = _lazy_import_predgpi_main()

    # 构建参数列表|Build argument list
    args = ["predgpi.py"]
    args.extend(["-i", input])
    args.extend(["-o", output_dir])
    if predgpi_home != _DEFAULT_PREDGPI_HOME:
        args.extend(["--predgpi-home", predgpi_home])
    if conservative:
        args.append("--conservative")
    if prefix:
        args.extend(["--prefix", prefix])

    original_argv = sys.argv
    sys.argv = args
    try:
        predgpi_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
