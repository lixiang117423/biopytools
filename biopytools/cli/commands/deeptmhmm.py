"""
DeepTMHMM命令|DeepTMHMM Command
"""

import click
import sys
import os


def _lazy_import_deeptmhmm_main():
    """延迟加载deeptmhmm主函数|Lazy load deeptmhmm main function"""
    try:
        from ...deeptmhmm.main import main as deeptmhmm_main
        return deeptmhmm_main
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
    short_help='DeepTMHMM跨膜螺旋/信号肽预测|DeepTMHMM TM helix & signal peptide prediction',
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
@click.option('--prefix',
              default=None,
              help='输出文件前缀(默认输入文件名)|Output file prefix')
@click.option('--conda-env',
              default='deeptmhmm_v.1.0',
              help='conda环境名|conda env name')
@click.option('--deeptmhmm-dir',
              default='~/software/deeptmhmm/DeepTMHMM-Academic-License-v1.0',
              help='DeepTMHMM安装目录|DeepTMHMM install directory')
def deeptmhmm(input, output_dir, prefix, conda_env, deeptmhmm_dir):
    """
    DeepTMHMM跨膜螺旋/信号肽预测|DeepTMHMM TM helix & signal peptide prediction

    示例|Examples: biopytools deeptmhmm -i proteins.fa -o output_dir
    """
    deeptmhmm_main = _lazy_import_deeptmhmm_main()

    args = ['deeptmhmm.py']
    args.extend(['-i', input])
    args.extend(['-o', output_dir])

    if prefix:
        args.extend(['--prefix', prefix])
    if conda_env != 'deeptmhmm_v.1.0':
        args.extend(['--conda-env', conda_env])
    if deeptmhmm_dir != '~/software/deeptmhmm/DeepTMHMM-Academic-License-v1.0':
        args.extend(['--deeptmhmm-dir', deeptmhmm_dir])

    original_argv = sys.argv
    sys.argv = args

    try:
        deeptmhmm_main()
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
