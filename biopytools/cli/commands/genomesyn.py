"""
🧠 基因组共线性分析命令 | Genome Synteny Analysis Command
"""

import click
import sys
# In your actual project structure, you would use a relative import like this:
# from ...genomesyn.main import main as genome_syn_main

# --- For this snippet to be self-contained, we'll include a placeholder main ---
# --- In your project, delete this placeholder and use the import above ---
# START: Placeholder for demonstration
def get_original_main_for_demo():
    def main_placeholder():
        print("--- 🚀 Original main function called (simulated) ---")
        print(f"Received sys.argv: {sys.argv}")
        # argparse would parse sys.argv here and run the analysis
        print("--- ✅ Analysis finished (simulated) ---")
    return main_placeholder
genome_syn_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--sample-map', '-s',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='📂 样本映射文件 (genome_file\\tgenome_name) | Sample mapping file.')
@click.option('--config', '-c',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='📋 配置文件 (.xlsx or .yaml) | Configuration file.')
@click.option('--output-dir', '-o',
              default='./genome_syn_output',
              show_default=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📁 输出目录 | Output directory.')
@click.option('--generate-config',
              is_flag=True,
              help='📝 仅生成配置文件 | Generate configuration file only.')
@click.option('--aligner', '-a',
              type=click.Choice(['minimap2', 'mcscanx', 'syri', 'mummer'], case_sensitive=False),
              default='minimap2',
              show_default=True,
              help='🔧 比对器类型 | Aligner type.')
@click.option('--alignment-mode',
              type=click.Choice(['chain', 'star', 'all_vs_all'], case_sensitive=False),
              default='chain',
              show_default=True,
              help='🔗 比对模式 | Alignment mode.')
@click.option('--threads', '-t',
              type=int,
              default=16,
              show_default=True,
              help='⚡ 线程数 | Number of threads.')
@click.option('--min-length',
              type=int,
              default=5000,
              show_default=True,
              help='📏 最小比对长度 | Minimum alignment length.')
@click.option('--chromosome',
              type=str,
              help='🧬 指定分析的染色体 (e.g., "1,2,3" or "1-5") | Specify chromosomes to analyze.')
@click.option('--canvas-width',
              type=int,
              help='📐 画布宽度 (自动计算) | Canvas width (auto-calculated).')
@click.option('--canvas-height',
              type=int,
              help='📐 画布高度 (自动计算) | Canvas height (auto-calculated).')
@click.option('--output-formats',
              multiple=True,
              type=click.Choice(['svg', 'png'], case_sensitive=False),
              default=['svg', 'png'],
              show_default=True,
              help='📄 输出格式 (可多次使用) | Output formats (can be used multiple times).')
def genomesyn(sample_map, config, output_dir, generate_config, aligner,
              alignment_mode, threads, min_length, chromosome,
              canvas_width, canvas_height, output_formats):
    """
    基因组共线性可视化工具.

    一个强大的工具，用于执行、配置和可视化多个基因组之间的共线性关系。
    
    🌟 示例 | Examples:
    
    \b
    # 📋 1. 从样本表生成配置文件
    biopytools genomesyn --sample-map genomes.tsv --generate-config -o ./output
    
    \b
    # 🚀 2. 使用生成的配置文件运行分析
    biopytools genomesyn --config ./output/genomes_config.xlsx -o ./output
    
    \b
    # ⚡ 3. (可选) 从样本表一步完成分析
    biopytools genomesyn --sample-map genomes.tsv -o ./output -t 32
    
    \b
    # 🧬 4. 仅分析特定染色体
    biopytools genomesyn -c config.xlsx -o ./output --chromosome "1,2,3"
    """
    
    # 验证互斥参数 | Validate mutually exclusive arguments
    if not sample_map and not config:
        raise click.UsageError("❌ 必须提供 '--sample-map' 或 '--config' 中的一个 | Either '--sample-map' or '--config' must be provided.")
    if sample_map and config:
        raise click.UsageError("❌ '--sample-map' 和 '--config' 不能同时使用 | '--sample-map' and '--config' are mutually exclusive.")

    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'genomesyn']

    # 必需的互斥参数 | Required mutually exclusive parameters
    if sample_map:
        args.extend(['--sample-map', sample_map])
    if config:
        args.extend(['--config', config])
        
    # 可选参数 ⚙️ | Optional parameters
    if output_dir != './genome_syn_output':
        args.extend(['-o', output_dir])
    else: # 总是传递输出目录，即使是默认值
        args.extend(['-o', output_dir])
        
    if generate_config:
        args.append('--generate-config')
        
    if aligner != 'minimap2':
        args.extend(['-a', aligner])
        
    if alignment_mode != 'chain':
        args.extend(['--alignment-mode', alignment_mode])
        
    if threads != 16:
        args.extend(['-t', str(threads)])
        
    if min_length != 5000:
        args.extend(['--min-length', str(min_length)])
        
    if chromosome:
        args.extend(['--chromosome', chromosome])
        
    if canvas_width:
        args.extend(['--canvas-width', str(canvas_width)])
        
    if canvas_height:
        args.extend(['--canvas-height', str(canvas_height)])
        
    # 处理多值参数 | Handle multiple value parameter
    # 只有当它与默认值不同时才添加
    if sorted(output_formats) != sorted(['svg', 'png']):
        # argparse的nargs='+'期望 `--flag val1 val2`
        args.append('--output-formats')
        args.extend(output_formats)

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        genome_syn_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        if e.code != 0:
            click.secho(f"❌ 脚本执行被终止，退出码: {e.code}", fg='red', err=True)
        sys.exit(e.code)
    except Exception as e:
        click.secho(f"💥 发生未知错误 | An unexpected error occurred: {e}", fg='red', err=True)
        sys.exit(1)
    finally:
        # 无论如何都要恢复原始的 sys.argv | Restore original sys.argv regardless of outcome
        sys.argv = original_argv

# 如果直接运行此文件用于测试 | If running this file directly for testing
if __name__ == '__main__':
    # 模拟命令行调用，例如: python your_script.py -s your_map.tsv -o ./out
    genomesyn()