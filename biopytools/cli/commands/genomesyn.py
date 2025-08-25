"""
基因组共线性分析命令 | Genome Synteny Analysis Command
优化版本：使用懒加载解决响应速度问题
"""

import click
import sys
import os


def _lazy_import_genome_syn_main():
    """懒加载genome_syn main函数 | Lazy load genome_syn main function"""
    try:
        from ...genome_syn.main import main as genome_syn_main
        return genome_syn_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在 | File does not exist: {file_path}")
    return file_path


@click.command(short_help='基因组共线性分析',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--sample-map', '-s',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📂 样本映射文件 | Sample mapping file (tab-separated: genome_file\\tgenome_name)')
@click.option('--config', '-c',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📋 配置文件 | Configuration file (.xlsx or .yaml)')
@click.option('--output-dir', '-o',
              default='./genome_syn_output',
              type=click.Path(),
              help='📁 输出目录 | Output directory (default: ./genome_syn_output)')
@click.option('--generate-config',
              is_flag=True,
              help='📋 仅生成配置文件 | Generate configuration file only')
@click.option('--aligner', '-a',
              default='minimap2',
              type=click.Choice(['minimap2', 'mcscanx', 'syri', 'mummer']),
              help='🔧 比对器类型 | Aligner type (default: minimap2)')
@click.option('--alignment-mode',
              default='chain',
              type=click.Choice(['chain', 'star', 'all_vs_all']),
              help='🔗 比对模式 | Alignment mode (default: chain)')
@click.option('--threads', '-t',
              default=16,
              type=int,
              help='⚡ 线程数 | Number of threads (default: 16)')
@click.option('--min-length',
              default=5000,
              type=int,
              help='📏 最小比对长度 | Minimum alignment length (default: 5000)')
@click.option('--chromosome',
              type=str,
              help='🧬 指定要分析的染色体 | Specify chromosomes to analyze (e.g., "1,2,3" or "1-5" or "1")')
@click.option('--canvas-width',
              type=int,
              help='📏 画布宽度 | Canvas width (auto-calculated if not specified)')
@click.option('--canvas-height',
              type=int,
              help='📏 画布高度 | Canvas height (auto-calculated if not specified)')
@click.option('--output-formats',
              multiple=True,
              default=['svg', 'png'],
              type=click.Choice(['svg', 'png']),
              help='📄 输出格式 | Output formats (default: svg png)')
def genomesyn(sample_map, config, output_dir, generate_config, aligner, 
               alignment_mode, threads, min_length, chromosome, canvas_width, 
               canvas_height, output_formats):
    """
    基因组共线性可视化工具
    
    从基因组序列文件进行共线性分析并生成可视化图形。
    支持多种比对器和可视化格式，可灵活配置分析参数。
    
    示例 | Examples:
    
    \b
    # 📋 生成配置文件
    biopytools genomesyn --sample-map genomes.tsv --generate-config --output-dir ./output
    
    \b
    # 🚀 运行分析
    biopytools genomesyn --config genomes_config.xlsx --output-dir ./output
    
    \b
    # ⚡ 一步完成
    biopytools genomesyn --sample-map genomes.tsv --output-dir ./output
    
    \b
    # 🧬 分析特定染色体
    biopytools genomesyn --sample-map genomes.tsv --output-dir ./output --chromosome "1,2,3"
    
    \b
    # 📊 分析染色体范围
    biopytools genomesyn --sample-map genomes.tsv --output-dir ./output --chromosome "1-5"
    
    \b
    # 🎨 自定义可视化
    biopytools genomesyn --sample-map genomes.tsv --output-dir ./output \\
        --canvas-width 1200 --canvas-height 800 --output-formats svg png
    
    \b
    # ⚙️ 高级设置
    biopytools genomesyn --sample-map genomes.tsv --output-dir ./output \\
        --aligner minimap2 --alignment-mode chain --threads 32 --min-length 10000
    """
    
    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    genome_syn_main = _lazy_import_genome_syn_main()
    
    # 验证输入参数的互斥性
    if not sample_map and not config:
        raise click.ClickException("❌ 必须指定 --sample-map 或 --config 参数之一 | Must specify either --sample-map or --config")
    
    if sample_map and config:
        raise click.ClickException("❌ --sample-map 和 --config 参数不能同时使用 | Cannot use both --sample-map and --config")
    
    # 构建参数列表传递给原始main函数
    args = ['genomesyn.py']
    
    # 输入文件参数（互斥）
    if sample_map:
        args.extend(['--sample-map', sample_map])
    if config:
        args.extend(['--config', config])
    
    # 输出参数
    args.extend(['--output-dir', output_dir])
    
    # 模式参数
    if generate_config:
        args.append('--generate-config')
    
    # 比对参数
    if aligner != 'minimap2':
        args.extend(['--aligner', aligner])
    
    if alignment_mode != 'chain':
        args.extend(['--alignment-mode', alignment_mode])
    
    if threads != 16:
        args.extend(['--threads', str(threads)])
    
    if min_length != 5000:
        args.extend(['--min-length', str(min_length)])
    
    # 染色体过滤参数
    if chromosome:
        args.extend(['--chromosome', chromosome])
    
    # 可视化参数
    if canvas_width:
        args.extend(['--canvas-width', str(canvas_width)])
    
    if canvas_height:
        args.extend(['--canvas-height', str(canvas_height)])
    
    # 输出格式参数（只在非默认值时添加）
    if output_formats and set(output_formats) != {'svg', 'png'}:
        args.extend(['--output-formats'] + list(output_formats))
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        genome_syn_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n⚠️ 用户中断操作 | User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 运行错误 | Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv