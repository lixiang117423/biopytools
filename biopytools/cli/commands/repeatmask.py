"""
重复序列屏蔽命令|Repeat Masking Command
"""

import click
import sys
import os


def _lazy_import_repeatmask_main():
    """延迟加载repeatmask主函数|Lazy load repeatmask main function"""
    try:
        from ...repeatmask.main import main as repeatmask_main
        return repeatmask_main
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
    short_help='重复序列屏蔽工具|Repeat masking tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--genome',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              help='基因组FASTA文件|Genome FASTA file')
@click.option('-o', '--output-dir',
              default='./repeatmask_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('-s', '--species',
              default=None,
              help='物种名称|Species name for Dfam/Repbase database')
@click.option('--repeatmodeler-path',
              default='~/miniforge3/envs/repeatmodeler_v.2.0.7/bin/RepeatModeler',
              show_default=True,
              help='RepeatModeler路径|RepeatModeler path')
@click.option('--repeatmasker-path',
              default='~/miniforge3/envs/repeat_identiy/bin/RepeatMasker',
              show_default=True,
              help='RepeatMasker路径|RepeatMasker path')
@click.option('--builddatabase-path',
              default='BuildDatabase',
              show_default=True,
              help='BuildDatabase路径|BuildDatabase path')
@click.option('--skip-modeler',
              is_flag=True,
              help='跳过RepeatModeler|Skip RepeatModeler step')
@click.option('--no-dfam',
              is_flag=True,
              help='不使用Dfam数据库|Do not use Dfam database')
@click.option('--masking-mode',
              type=click.Choice(['soft', 'hard', 'x']),
              default='soft',
              show_default=True,
              help='屏蔽模式|Masking mode')
@click.option('--no-ltr',
              is_flag=True,
              help='不运行LTR结构发现|Do not run LTR structural discovery')
@click.option('--modeler-quick',
              is_flag=True,
              help='RepeatModeler快速模式|RepeatModeler quick mode')
def repeatmask(genome, output_dir, threads, species,
               repeatmodeler_path, repeatmasker_path, builddatabase_path,
               skip_modeler, no_dfam, masking_mode, no_ltr, modeler_quick):
    """
    重复序列屏蔽工具|Repeat Masking Tool

    使用RepeatModeler和RepeatMasker进行基因组重复序列识别和屏蔽|Use RepeatModeler and RepeatMasker for genome repeat identification and masking

    示例|Examples: biopytools repeatmask -i genome.fa -o output
    """

    # 延迟加载|Lazy loading
    repeatmask_main = _lazy_import_repeatmask_main()

    # 构建参数列表|Build argument list
    args = ['repeatmask.py']

    # 必需参数|Required parameters
    args.extend(['-i', genome])

    # 可选参数|Optional parameters
    if output_dir != './repeatmask_output':
        args.extend(['-o', output_dir])

    if threads != 12:
        args.extend(['-t', str(threads)])

    if species:
        args.extend(['-s', species])

    if repeatmodeler_path != '~/miniforge3/envs/repeatmodeler_v.2.0.7/bin/RepeatModeler':
        args.extend(['--repeatmodeler-path', repeatmodeler_path])

    if repeatmasker_path != '~/miniforge3/envs/repeat_identiy/bin/RepeatMasker':
        args.extend(['--repeatmasker-path', repeatmasker_path])

    if builddatabase_path != 'BuildDatabase':
        args.extend(['--builddatabase-path', builddatabase_path])

    # 布尔选项|Boolean options
    if skip_modeler:
        args.append('--skip-modeler')

    if no_dfam:
        args.append('--no-dfam')

    if masking_mode != 'soft':
        args.extend(['--masking-mode', masking_mode])

    if no_ltr:
        args.append('--no-ltr')

    if modeler_quick:
        args.append('--modeler-quick')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        repeatmask_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
