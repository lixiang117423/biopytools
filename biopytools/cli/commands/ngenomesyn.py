"""
基因组共线性分析|Genome Synteny Analysis Command
"""

import click
import sys
import os


def _lazy_import_ngenomesyn_main():
    """延迟加载ngenomesyn主函数|Lazy load ngenomesyn main function"""
    try:
        from ...ngenomesyn.main import main as ngenomesyn_main
        return ngenomesyn_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='基因组共线性分析工具|Genome synteny analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--sample-map', '-s',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='样本映射文件|Sample mapping file (genome_file\\tgenome_name)')
@click.option('--config', '-c',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='配置文件|Configuration file')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--aligner', '-a',
              default='minimap2',
              show_default=True,
              type=click.Choice(['minimap2', 'mummer']),
              help='比对器类型|Aligner type')
@click.option('--threads', '-t',
              default=12,
              show_default=True,
              type=int,
              help='线程数|Number of threads')
@click.option('--min-length',
              default=5000,
              show_default=True,
              type=int,
              help='最小比对长度|Minimum alignment length')
@click.option('--minimap-preset',
              default='asm5',
              show_default=True,
              help='Minimap2预设模式|Minimap2 preset')
@click.option('--mummer-match-type',
              default='mumreference',
              show_default=True,
              type=click.Choice(['mum', 'mumreference', 'maxmatch']),
              help='MUMmer匹配类型|MUMmer match type')
@click.option('--mummer-min-match',
              default=20,
              show_default=True,
              type=int,
              help='MUMmer最小匹配长度|MUMmer min match')
@click.option('--chromosome',
              type=str,
              help='指定染色体|Specify chromosomes (e.g., "1,2,3" or "1-5")')
@click.option('--output-formats',
              multiple=True,
              default=['svg', 'png'],
              show_default=True,
              type=click.Choice(['svg', 'png']),
              help='输出格式|Output formats')
@click.option('--ngenomesyn-bin',
              help='NGenomeSyn二进制文件路径|NGenomeSyn binary path')
@click.option('--use-syri',
              is_flag=True,
              default=False,
              help='使用SyRI进行结构变异分析|Use SyRI for structural variation analysis')
@click.option('--syri-bin',
              help='SyRI二进制文件路径|SyRI binary path')
def ngenomesyn(sample_map, config, output, aligner, threads, min_length,
               minimap_preset, mummer_match_type, mummer_min_match,
               chromosome, output_formats, ngenomesyn_bin, use_syri, syri_bin):
    """
    基因组共线性分析工具|Genome Synteny Analysis Tool

    基于NGenomeSyn进行基因组共线性分析和可视化
    Perform genome synteny analysis and visualization using NGenomeSyn

    示例|Examples: biopytools ngenomesyn -s samples.txt -o output_dir
    """

    # 延迟加载|Lazy load
    ngenomesyn_main = _lazy_import_ngenomesyn_main()

    # 验证参数|Validate parameters
    if not sample_map and not config:
        raise click.ClickException("必须指定--sample-map或--config|Must specify either --sample-map or --config")

    if sample_map and config:
        raise click.ClickException("不能同时使用--sample-map和--config|Cannot use both --sample-map and --config")

    # 构建主函数参数列表|Build argument list for main function
    args = ['ngenomesyn.py']

    # 输入参数|Input parameters
    if sample_map:
        args.extend(['-s', sample_map])
    if config:
        args.extend(['-c', config])

    # 输出参数|Output parameters
    args.extend(['-o', output])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if aligner != 'minimap2':
        args.extend(['--aligner', aligner])

    if threads != 16:
        args.extend(['--threads', str(threads)])

    if min_length != 5000:
        args.extend(['--min-length', str(min_length)])

    if minimap_preset != 'asm5':
        args.extend(['--minimap-preset', minimap_preset])

    if mummer_match_type != 'mumreference':
        args.extend(['--mummer-match-type', mummer_match_type])

    if mummer_min_match != 20:
        args.extend(['--mummer-min-match', str(mummer_min_match)])

    # 染色体参数|Chromosome parameter
    if chromosome:
        args.extend(['--chromosome', chromosome])

    # 输出格式|Output formats
    if output_formats and set(output_formats) != {'svg', 'png'}:
        args.extend(['--output-formats'] + list(output_formats))

    if ngenomesyn_bin:
        args.extend(['--ngenomesyn-bin', ngenomesyn_bin])

    # SyRI参数|SyRI parameters
    if use_syri:
        args.append('--use-syri')

    if syri_bin:
        args.extend(['--syri-bin', syri_bin])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        ngenomesyn_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
