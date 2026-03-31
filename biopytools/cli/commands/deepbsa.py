"""
DeepBSA批量分析命令|DeepBSA Batch Analysis Command
"""

import click
import sys
import os


def _lazy_import_deepbsa_main():
    """延迟加载deepbsa主函数|Lazy load deepbsa main function"""
    try:
        from ...deepbsa.main import main as deepbsa_main
        return deepbsa_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if '-h' not in sys.argv and '--help' not in sys.argv:
        if file_path and not os.path.exists(file_path):
            raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='DeepBSA批量分析工具|DeepBSA batch analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入文件路径（VCF/CSV格式）|Input file path (VCF/CSV format)')
@click.option('-m', '--methods',
              default=None,
              help='要运行的方法，逗号分隔（默认：全部）|Methods to run, comma-separated (default: all)\n可用方法|Available: DL,K,ED4,SNP,SmoothG,SmoothLOD,Ridit')
@click.option('-o', '--output-dir',
              default='deepbsa_results',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('-p', '--parallel',
              is_flag=True,
              default=True,
              help='并行运行所有方法（默认）|Run all methods in parallel (default)')
@click.option('--no-parallel',
              is_flag=True,
              help='串行运行所有方法|Run all methods serially')
@click.option('-n', '--no-auto-clean',
              is_flag=True,
              help='不自动清理VCF注释行|Don\'t auto clean VCF comment lines')
@click.option('-k', '--keep-clean',
              is_flag=True,
              help='保留清理后的文件|Keep cleaned file')
@click.option('--deepbsa-path',
              default='~/software/DeepBSA/DeepBSA_linux_v1.4/bin/main.py',
              show_default=True,
              help='DeepBSA主程序路径|DeepBSA main script path')
@click.option('--conda-env',
              default='/share/org/YZWL/yzwl_lixg/miniforge3/envs/DeepBSA',
              show_default=True,
              help='Conda环境路径|Conda environment path')
@click.option('-v', '--verbose',
              is_flag=True,
              help='详细输出|Verbose output')
def deepbsa(input_file, methods, output_dir, parallel, no_parallel, no_auto_clean, keep_clean,
            deepbsa_path, conda_env, verbose):
    """
    DeepBSA批量分析工具|DeepBSA Batch Analysis Tool

    批量运行多种BSA分析方法（DL、K、ED4、SNP、SmoothG、SmoothLOD、Ridit）|Run multiple BSA analysis methods (DL, K, ED4, SNP, SmoothG, SmoothLOD, Ridit)

    示例|Example: biopytools deepbsa -i filtered.vcf -m DL,K,ED4

    默认并行运行多个方法，使用--no-parallel串行运行|Run multiple methods in parallel by default, use --no-parallel for serial execution

    可用方法|Available methods: DL, K, ED4, SNP, SmoothG, SmoothLOD, Ridit
    """

    # 延迟加载|Lazy loading
    deepbsa_main = _lazy_import_deepbsa_main()

    # 构建参数列表|Build argument list
    args = ['deepbsa.py']

    # 必需参数|Required parameters
    args.extend(['-i', input_file])

    # 方法选择|Method selection
    if methods:
        args.extend(['-m', methods])

    # 输出参数|Output parameters
    if output_dir != 'deepbsa_results':
        args.extend(['-o', output_dir])

    # 运行模式|Running mode
    # 处理parallel和no_parallel的逻辑
    actual_parallel = parallel and not no_parallel

    if actual_parallel:
        args.append('-p')

    # 预处理参数|Preprocessing parameters
    if no_auto_clean:
        args.append('-n')

    if keep_clean:
        args.append('-k')

    # 工具路径|Tool paths
    if deepbsa_path != '~/software/DeepBSA/DeepBSA_linux_v1.4/bin/main.py':
        args.extend(['--deepbsa-path', deepbsa_path])

    if conda_env != '/share/org/YZWL/yzwl_lixg/miniforge3/envs/DeepBSA':
        args.extend(['--conda-env', conda_env])

    # 其他参数|Other parameters
    if verbose:
        args.append('-v')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        deepbsa_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


