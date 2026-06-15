"""
基因组共线性分析命令|Genome Synteny Analysis Command
"""

import click
import sys
import os


def _lazy_import_genome_syn_main():
    """延迟加载genome_syn主函数|Lazy load genome_syn main function"""
    try:
        from ...genomesyn.main import main as genome_syn_main
        return genome_syn_main
    except ImportError as e:
        click.echo(f"Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在性(仅在非帮助模式下)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='基因组共线性分析|Genome Synteny Analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--ref',
              type=click.Path(exists=True),
              help='参考基因组|Reference genome file')
@click.option('-I', '--query',
              multiple=True,
              type=click.Path(exists=True),
              help='比对基因组（可重复）|Query genome file(s), can be repeated')
@click.option('--sample-map', '-s',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='样本映射文件|Sample mapping file')
@click.option('--config', '-c',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='配置文件(.xlsx或.yaml)|Configuration file (.xlsx or .yaml)')
@click.option('--output-dir', '-o',
              default='./genome_syn_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--generate-config',
              is_flag=True,
              help='仅生成配置文件|Generate configuration file only')
@click.option('--aligner', '-a',
              default='minimap2',
              show_default=True,
              type=click.Choice(['minimap2', 'mcscanx', 'syri', 'mummer']),
              help='比对工具|Aligner type')
@click.option('--alignment-mode',
              default='chain',
              show_default=True,
              type=click.Choice(['chain', 'star', 'all_vs_all']),
              help='比对模式|Alignment mode')
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
@click.option('--no-minimap2-preset',
              is_flag=True,
              default=False,
              help='不使用minimap2的-x预设参数|Disable minimap2 -x preset parameter')
@click.option('--chromosome',
              type=str,
              help='指定染色体|Specify chromosomes to analyze')
@click.option('--canvas-width',
              type=int,
              help='画布宽度|Canvas width')
@click.option('--canvas-height',
              type=int,
              help='画布高度|Canvas height')
@click.option('--output-formats',
              multiple=True,
              default=['svg', 'png'],
              show_default=True,
              type=click.Choice(['svg', 'png']),
              help='输出格式|Output formats')
@click.option('--regions',
              type=click.Path(exists=True),
              help='区域标注文件|Special region annotation file')
def genomesyn(ref, query, sample_map, config, output_dir, generate_config, aligner,
               alignment_mode, threads, min_length, no_minimap2_preset, chromosome,
               canvas_width, canvas_height, output_formats, regions):
    """
    基因组共线性分析工具|Genome Synteny Analysis Tool

    多基因组共线性可视化和分析|Multi-genome synteny visualization and analysis

    示例|Examples: biopytools genomesyn -i ref.fa -I query.fa -o output
    """

    # 延迟加载|Lazy loading: import only when actually called
    genome_syn_main = _lazy_import_genome_syn_main()

    # 验证参数|Validate parameters
    if config and (ref or query or sample_map):
        raise click.ClickException("--config不能与其他参数同时使用|--config cannot be used with other parameters")

    if not config and not sample_map and not (ref and query):
        raise click.ClickException("必须指定--config或--sample-map或同时指定-i和-I|Must specify --config or --sample-map or both -i and -I")

    if sample_map and (ref or query):
        raise click.ClickException("不能同时使用--sample-map和-i/-I|Cannot use both --sample-map and -i/-I")

    # 如果使用 -i/-I 模式，生成临时 sample_map 文件
    temp_sample_map = None
    if ref and query:
        import os

        # 合并所有基因组：ref + queries
        all_genomes = [ref] + list(query)

        # 生成临时 sample_map 文件
        temp_sample_map = os.path.join(output_dir, '.temp_sample_map.txt')

        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)

        # 从文件名推断基因组名
        def get_genome_name(file_path):
            return os.path.splitext(os.path.basename(file_path))[0]

        # 写入临时文件
        with open(temp_sample_map, 'w') as f:
            for genome_file in all_genomes:
                genome_name = get_genome_name(genome_file)
                f.write(f"{genome_file}\t{genome_name}\n")

    # 构建参数列表|Build argument list
    args = ['genomesyn.py']

    # 样本映射或配置文件|Sample map or config file
    if sample_map:
        args.extend(['--sample-map', sample_map])
    elif temp_sample_map:
        args.extend(['--sample-map', temp_sample_map])
    if config:
        args.extend(['--config', config])

    # 输出目录|Output directory
    args.extend(['--output-dir', output_dir])

    # 配置生成|Config generation
    if generate_config:
        args.append('--generate-config')

    # 比对参数|Alignment parameters
    if aligner != 'minimap2':
        args.extend(['--aligner', aligner])

    if alignment_mode != 'chain':
        args.extend(['--alignment-mode', alignment_mode])

    if threads != 16:
        args.extend(['--threads', str(threads)])

    if min_length != 5000:
        args.extend(['--min-length', str(min_length)])

    if no_minimap2_preset:
        args.append('--no-minimap2-preset')

    # 染色体过滤|Chromosome filter
    if chromosome:
        args.extend(['--chromosome', chromosome])

    # 可视化参数|Visualization parameters
    if canvas_width:
        args.extend(['--canvas-width', str(canvas_width)])

    if canvas_height:
        args.extend(['--canvas-height', str(canvas_height)])

    # 区域标注|Region annotation
    if regions:
        args.extend(['--regions', regions])

    # 输出格式|Output formats
    if output_formats and set(output_formats) != {'svg', 'png'}:
        args.extend(['--output-formats'] + list(output_formats))

    # 保存和恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        genome_syn_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"执行失败|Execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
