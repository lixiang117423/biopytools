"""
TMHMM命令|TMHMM Command
"""

import click
import sys
import os


def _lazy_import_tmhmm_main():
    """延迟加载tmhmm主函数|Lazy load tmhmm main function"""
    try:
        from ...tmhmm.main import main as tmhmm_main
        return tmhmm_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在|Validate path exists"""
    if _is_help_request():
        return path
    if not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='TMHMM跨膜螺旋预测|TMHMM transmembrane helix prediction',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入蛋白质FASTA文件|Input protein FASTA file')
@click.option('-o', '--output-dir',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--plot',
              is_flag=True,
              default=False,
              help='生成图形(默认-noplot)|Generate plots (default -noplot)')
@click.option('--prefix',
              default=None,
              help='输出文件前缀|Output file prefix')
def tmhmm(input, output_dir, plot, prefix):
    """
    TMHMM跨膜螺旋预测|TMHMM transmembrane helix prediction

    示例|Examples: biopytools tmhmm -i proteins.fa -o output_dir
    """
    tmhmm_main = _lazy_import_tmhmm_main()

    args = ['tmhmm.py']
    args.extend(['-i', input])
    args.extend(['-o', output_dir])

    if plot:
        args.append('--plot')
    if prefix:
        args.extend(['--prefix', prefix])

    original_argv = sys.argv
    sys.argv = args

    try:
        tmhmm_main()
    except SystemExit as e:
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
