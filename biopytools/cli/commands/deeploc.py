"""
DeepLoc 2.1蛋白质亚细胞定位预测命令|DeepLoc 2.1 Protein Subcellular Localization Prediction Command
"""

import click
import sys
import os


def _lazy_import_deeploc_main():
    """延迟加载deeploc主函数|Lazy load deeploc main function"""
    try:
        from ...deeploc.main import main as deeploc_main
        return deeploc_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='DeepLoc 2.1蛋白质亚细胞定位预测工具|DeepLoc 2.1 protein subcellular localization prediction tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTA文件(蛋白质序列)|Input FASTA file (protein sequences)')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--model',
              default='Fast',
              type=click.Choice(['Fast', 'Accurate'], case_sensitive=False),
              show_default=True,
              help='预测模型|Prediction model')
@click.option('--device',
              default='cpu',
              type=click.Choice(['cpu', 'cuda', 'mps'], case_sensitive=False),
              show_default=True,
              help='计算设备|Compute device')
@click.option('--singularity-image',
              default='~/software/singularity/deeploc2.1_latest.sif',
              show_default=True,
              help='Singularity镜像路径|Singularity image path')
@click.option('--database-dir',
              default='~/software/deeploc/database',
              show_default=True,
              type=click.Path(),
              help='数据库目录|Database directory')
@click.option('--singularity-exec',
              default='~/miniforge3/envs/singularity_v.3.8.7/bin/singularity',
              show_default=True,
              help='Singularity可执行文件路径|Singularity executable path')
@click.option('--plot',
              is_flag=True,
              help='绘制attention图|Plot attention values')
def deeploc(input, output_dir, model, device, singularity_image, database_dir,
            singularity_exec, plot):
    """
    DeepLoc 2.1蛋白质亚细胞定位预测工具|DeepLoc 2.1 Protein Subcellular Localization Prediction Tool

    预测蛋白质的亚细胞定位和膜关联性|Predict protein subcellular localization and membrane association

    示例|Examples: biopytools deeploc -i proteins.faa -o output_dir
    """

    # 延迟加载|Lazy loading
    deeploc_main = _lazy_import_deeploc_main()

    # 构建参数列表|Build argument list
    args = ['deeploc.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output_dir])

    # 可选参数|Optional parameters
    if model != 'Fast':
        args.extend(['--model', model])

    if device != 'cpu':
        args.extend(['--device', device])

    if singularity_image != '~/software/singularity/deeploc2.1_latest.sif':
        args.extend(['--singularity-image', singularity_image])

    if str(database_dir) != '~/software/deeploc/database':
        args.extend(['--database-dir', database_dir])

    if str(singularity_exec) != '~/miniforge3/envs/singularity_v.3.8.7/bin/singularity':
        args.extend(['--singularity-exec', singularity_exec])

    if plot:
        args.append('--plot')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        deeploc_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
