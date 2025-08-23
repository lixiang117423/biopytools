"""
🧬 HiFiasm基因组组装命令 | HiFiasm Genome Assembly Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...hifiasm.main import main as hifiasm_main

# --- Placeholder for the original main function to make this snippet runnable ---
# --- In your project, delete this section and use the import statement above ---
# START: Placeholder for demonstration
def get_original_main_for_demo():
    def main_placeholder():
        print("--- 🚀 Original main function called (simulated) ---")
        print(f"Received sys.argv: {sys.argv}")
        # The original main() would create a parser and run the analysis
        print("--- ✅ Analysis finished (simulated) ---")
    return main_placeholder
hifiasm_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']), short_help = "HiFiasm基因组组装完整流水线")
# --- Required arguments ---
@click.option('--input-reads', '-i',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 输入HiFi测序数据文件 | Input HiFi sequencing data file.')
# --- Basic arguments ---
@click.option('--output-dir', '-o',
              default='./hifiasm_output', show_default=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📁 输出目录 | Output directory.')
@click.option('--prefix', '-p',
              default='sample', show_default=True,
              help='🏷️ 输出文件前缀 | Output file prefix.')
@click.option('--threads', '-t',
              type=int, default=32, show_default=True,
              help='🧵 线程数 | Number of threads.')
# --- HiFiasm assembly parameters ---
@click.option('--hg-size',
              default='auto', show_default=True,
              help='📏 基因组大小估计 (e.g., 1.4g) | Genome size estimation.')
@click.option('--purge-level', '-l',
              type=int, default=3, show_default=True,
              help='🧼 Purge级别 (0-3) | Purge level.')
@click.option('--purge-max',
              type=int, default=65, show_default=True,
              help='📈 最大purge覆盖度 | Maximum purge coverage.')
@click.option('--similarity-threshold', '-s',
              type=float, default=0.75, show_default=True,
              help='🎯 相似性阈值 | Similarity threshold.')
@click.option('--ont-reads',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 ONT长读长数据文件 | ONT long-read data file.')
@click.option('--hi-c-1',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🔗 Hi-C第一端数据文件 | Hi-C first-end data file.')
@click.option('--hi-c-2',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🔗 Hi-C第二端数据文件 | Hi-C second-end data file.')
@click.option('--extra-hifiasm-args',
              default='',
              help='⚙️ 额外的HiFiasm参数 (用引号括起来) | Additional HiFiasm arguments (in quotes).')
# --- Quality assessment parameters ---
@click.option('--skip-busco', is_flag=True, help='⏩ 跳过BUSCO质量评估 | Skip BUSCO quality assessment.')
@click.option('--busco-lineage', default='auto', show_default=True, help='🌳 BUSCO谱系数据集 | BUSCO lineage dataset.')
@click.option('--busco-mode', type=click.Choice(['genome', 'proteins', 'transcriptome']), default='genome', show_default=True, help='📊 BUSCO评估模式 | BUSCO assessment mode.')
@click.option('--skip-quast', is_flag=True, help='⏩ 跳过QUAST质量评估 | Skip QUAST quality assessment.')
@click.option('--reference-genome', type=click.Path(exists=True, dir_okay=False, resolve_path=True), help='🎯 参考基因组文件 (用于QUAST) | Reference genome file (for QUAST).')
# --- Analysis parameters ---
@click.option('--analyze-haplotypes', is_flag=True, help='🔬 分析单倍型差异 | Analyze haplotype differences.')
@click.option('--min-contig-length', type=int, default=1000, show_default=True, help='✂️ 最小contig长度过滤 | Minimum contig length filter.')
@click.option('--generate-plots', is_flag=True, help='📈 生成可视化图表 | Generate visualization plots.')
@click.option('--assembly-type', type=click.Choice(['auto', 'diploid', 'triploid', 'polyploid']), default='auto', show_default=True, help='🧬 组装类型 | Assembly type.')
# --- Output control parameters ---
@click.option('--keep-intermediate', is_flag=True, help='💾 保留中间文件 | Keep intermediate files.')
@click.option('--compress-output', is_flag=True, help='📦 压缩输出文件 | Compress output files.')
@click.option('--output-formats', multiple=True, type=click.Choice(['fasta', 'gfa', 'both']), default=['both'], show_default=True, help='📄 输出格式选择 | Output format selection.')
# --- System parameters ---
@click.option('--memory', type=int, default=64, show_default=True, help='🧠 内存大小(GB) | Memory size (GB).')
@click.option('--tmp-dir', default='/tmp', show_default=True, type=click.Path(file_okay=False), help='🗑️ 临时目录 | Temporary directory.')
@click.option('--max-runtime', type=int, default=48, show_default=True, help='⏰ 最大运行时间(小时) | Maximum runtime (hours).')
@click.option('--resume', is_flag=True, help='🔁 恢复中断的分析 | Resume interrupted analysis.')
# --- Tool paths parameters ---
@click.option('--hifiasm-path', default='hifiasm', show_default=True, help='🔧 HiFiasm软件路径 | HiFiasm software path.')
@click.option('--busco-path', default='busco', show_default=True, help='🔧 BUSCO软件路径 | BUSCO software path.')
@click.option('--quast-path', default='quast', show_default=True, help='🔧 QUAST软件路径 | QUAST software path.')
@click.option('--python-path', default='python3', show_default=True, help='🐍 Python解释器路径 | Python interpreter path.')
@click.option('--samtools-path', default='samtools', show_default=True, help='🔧 Samtools软件路径 | Samtools software path.')
# --- Database paths parameters ---
@click.option('--busco-db-path', type=click.Path(exists=True, file_okay=False), help='🗃️ BUSCO数据库路径 | BUSCO database path.')
@click.option('--busco-download-path', type=click.Path(file_okay=False), help='📥 BUSCO数据集下载路径 | BUSCO dataset download path.')
# --- Advanced parameters ---
@click.option('--debug', is_flag=True, help='🐞 启用调试模式 | Enable debug mode.')
@click.option('--verbose', '-v', count=True, help='🔊 详细输出模式 (-v, -vv) | Verbose output mode.')
@click.option('--log-level', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']), default='INFO', show_default=True, help='📜 日志级别 | Log level.')
@click.option('--config-file', type=click.Path(exists=True, dir_okay=False), help='🗒️ 配置文件路径 | Configuration file path.')
@click.option('--dry-run', is_flag=True, help='🧪 试运行模式 (不执行实际命令) | Dry run mode.')

def hifiasm(**kwargs):
    """
    HiFiasm基因组组装完整流水线.

    一个从HiFi reads到高质量、分单倍型基因组组装的端到端分析工具。
    
    示例 | Examples:
    
    \b
    # 🎯 基本用法 (HiFi reads -> 组装)
    biopytools hifiasm -i reads.fq.gz -o results -p sample1 -t 64
    
    \b
    # 🔗 使用Hi-C数据辅助分相
    biopytools hifiasm -i hifi.fq.gz -o results_hic -p sample2 \\
        --hi-c-1 hic_1.fq.gz --hi-c-2 hic_2.fq.gz
        
    \b
    # 🌳 运行BUSCO并生成图表
    biopytools hifiasm -i hifi.fq.gz -o results_busco -p sample3 \\
        --busco-lineage embryophyta_odb10 --generate-plots
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'hifiasm']
    
    # 遍历所有click传递的参数
    for key, value in kwargs.items():
        if value is None:
            continue
            
        param_name = '--' + key.replace('_', '-')
        
        # 处理布尔标志 | Handle boolean flags
        if isinstance(value, bool) and value:
            args.append(param_name)
        # 处理计数标志 | Handle count flags
        elif key == 'verbose' and value > 0:
            args.extend(['-v'] * value)
        # 处理多值参数 | Handle multiple value parameters
        elif isinstance(value, tuple) and value:
            # 检查是否与默认值不同
            default_val = hifiasm.params_by_name[key].default
            if sorted(value) != sorted(default_val):
                args.append(param_name)
                args.extend(value)
        # 处理非布尔、非默认值的常规参数 | Handle regular non-boolean, non-default parameters
        elif not isinstance(value, bool) and not isinstance(value, tuple):
            default_val = hifiasm.params_by_name[key].default
            # 只有当值不等于默认值时才添加，保持命令简洁
            # 对于必需参数（默认值为None），也添加
            if value != default_val or default_val is None:
                args.append(param_name)
                args.append(str(value))

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        hifiasm_main()
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
    hifiasm()