"""
🧬 基因组组装流程CLI包装器 | Genome Assembly Pipeline CLI Wrapper
高级优化版本：使用Click和懒加载以提高--help响应速度。
"""

import click
import sys
import os
from pathlib import Path

def _lazy_import_assembler_main():
    """懒加载基因组组装主函数 | Lazy load the genome assembler main function"""
    try:
        # 假设main.py在当前模块的父级目录
        # 如果结构是 genomeasm/main.py 和 genomeasm/cli.py, 使用 from .main import main
        from ...genomeasm import main as assembler_main
        return assembler_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: 无法加载组装主模块。请检查项目结构。", err=True)
        click.echo(f"   详细信息 | Details: {e}", err=True)
        sys.exit(1)

def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)

def _validate_input_dir(ctx, param, value):
    """验证输入目录是否存在（仅在非帮助模式下）| Validate input directory (only in non-help mode)"""
    if not _is_help_request() and not Path(value).is_dir():
        raise click.BadParameter(f"输入目录不存在或不是一个目录 | Input directory does not exist or is not a directory: {value}")
    return value

@click.command(
    name="assemble",
    short_help="🧬 运行基因组组装流程",
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120),
    help="""
    🧬 基因组组装工具 (多数据整合版本) | Genome Assembly Tool (Multi-data Integration Version)

    🌟 示例用法 | Examples:
    
    \b
    # 基础HiFi组装
    biopytools assemble -i raw_data/ -o assembly_results/
    
    \b
    # HiFi + Hi-C染色体级组装
    biopytools assemble -i data/ -o results/ --hic-strategy complete_juicer
    
    \b
    # 指定项目参数
    biopytools assemble -i input/ -o output/ -n my_genome --genome-size 3g -t 64
    
    \b
    # 使用简化Hi-C流程
    biopytools assemble -i data/ -o results/ --hic-strategy simplified_salsa2

    📚 支持的数据类型 | Supported Data Types:
      - HiFi: 高准确长读长数据 (必需) | High-accuracy long reads (required)
      - Hi-C: 染色体构象捕获数据 | Chromosome conformation capture data
      - ONT: Oxford Nanopore长读长数据 | Oxford Nanopore long reads
      - NGS: Illumina短读长数据 | Illumina short reads

    🔗 Hi-C处理策略 | Hi-C Processing Strategies:
      - complete_juicer: 完整Juicer + 3D-DNA流程 (最高质量)
      - standard_3ddna: 简化3D-DNA流程 (平衡质量与复杂度)
      - simplified_salsa2: SALSA2流程 (最简化)
    """
)
# --- 必需参数 | Required Arguments ---
@click.option('-i', '--input-dir', 
              required=True, 
              type=click.Path(exists=False, file_okay=False, dir_okay=True),
              callback=_validate_input_dir,
              help='输入数据目录 (自动检测文件类型)')

# --- 基本参数 | Basic Arguments ---
@click.option('-o', '--output-dir', default='./assembly_output', type=click.Path(), help='输出目录')
@click.option('-n', '--project-name', default='genome_assembly', help='项目名称')
@click.option('-t', '--threads', type=int, default=88, help='线程数')

# --- Hi-C参数 | Hi-C Arguments ---
@click.option('--hic-strategy', 
              type=click.Choice(['complete_juicer', 'standard_3ddna', 'simplified_salsa2']),
              default='complete_juicer', help='Hi-C处理策略')
@click.option('--restriction-enzyme', 
              type=click.Choice(['MboI', 'DpnII', 'HindIII', 'EcoRI']),
              default='MboI', help='限制性酶类型')
@click.option('--min-contig-size', type=int, default=15000, help='最小contig大小阈值')
@click.option('--edit-rounds', type=int, default=2, help='3D-DNA编辑轮数')

# --- 组装参数 | Assembly Arguments ---
@click.option('--genome-size', default='3g', help='预估基因组大小 (e.g., 3g, 500m)')
@click.option('--species-type', 
              type=click.Choice(['diploid', 'haploid', 'polyploid']),
              default='diploid', help='物种倍性')
@click.option('--telomere-motif', default='CCCTAA', help='端粒序列motif')
@click.option('--purge-level', type=click.Choice(['0', '1', '2', '3']), default='1', help='Purging级别')
@click.option('--purge-max', type=int, default=80, help='Purging覆盖度上限')
@click.option('--similarity-threshold', type=float, default=0.75, help='相似度阈值')
@click.option('--n-haplotypes', type=int, default=2, help='单倍型数量')

# --- 质量控制参数 | Quality Control Arguments ---
# @click.option('--skip-fastqc', default=True, help='跳过FastQC质量检查 (默认跳过，节省时间) | Skip FastQC quality check (default: skip to save time)')
@click.option('--skip-fastqc', default=True, help='跳过FastQC质量检查 (默认跳过，节省时间)')
@click.option('--min-hifi-coverage', type=int, default=30, help='最小HiFi覆盖度')
@click.option('--min-hic-coverage', type=int, default=50, help='最小Hi-C覆盖度')
@click.option('--min-mapping-rate', type=float, default=0.7, help='最小映射率')
@click.option('--busco-lineage', default='auto', help='BUSCO谱系数据库')

# --- 工具路径参数 | Tool Paths ---
@click.option('--hifiasm-path', default='hifiasm', help='Hifiasm程序路径')
@click.option('--bwa-path', default='bwa', help='BWA程序路径')
@click.option('--samtools-path', default='samtools', help='Samtools程序路径')
@click.option('--juicer-path', default='juicer.sh', help='Juicer脚本路径')
@click.option('--pipeline-3ddna', default='3d-dna/run-asm-pipeline.sh', help='3D-DNA pipeline路径')
@click.option('--juicer-tools', default='juicer_tools.jar', help='Juicer tools JAR路径')
@click.option('--salsa2-path', default='run_pipeline.py', help='SALSA2脚本路径')
def genomeasm(**kwargs):
    """
    🧬 基因组组装流程的Click包装器入口点。
    
    此函数负责收集所有命令行参数，并将它们转换为与原始
    `argparse` 脚本兼容的 `sys.argv` 列表，然后调用主逻辑。
    """
    # 🚀 懒加载：只有在实际调用时才导入主模块
    assembler_main = _lazy_import_assembler_main()
    
    # 🔄 构建参数列表以传递给原始main函数
    # 第一个参数是脚本名，与argparse的行为保持一致
    args = ['genomeasm.py']
    
    # 获取所有定义选项的默认值
    defaults = {param.name: param.default for param in genomeasm.params}
    
    for key, value in kwargs.items():
        # click将命令行选项的'-'替换为'_'
        cli_option = '--' + key.replace('_', '-')
        
        # 仅当用户提供的值不是默认值或该选项为必需选项时，才将其添加到参数列表
        # 'input_dir' 是必需的，所以它没有默认值 (None)
        if value is not None and value != defaults.get(key):
            args.extend([cli_option, str(value)])

    # 💾 保存并恢复sys.argv，确保隔离执行
    original_argv = sys.argv
    sys.argv = args
    
    # click.echo("🚀 启动基因组组装流程...")
    # click.echo(f"   执行命令: {' '.join(sys.argv)}")
    
    try:
        # 🚀 调用原始的main函数
        assembler_main()
    except SystemExit as e:
        # ✅ 处理程序正常退出 (e.g., sys.exit(0) or sys.exit(1))
        # click会自动处理，这里可以留空或记录日志
        if e.code != 0:
            click.secho(f"❌ 流程因错误而终止，退出代码: {e.code}", fg='red', err=True)
        sys.exit(e.code)
    except Exception as e:
        # 💥 捕获其他意外异常
        click.secho(f"💥 程序执行失败: {e}", fg='red', bold=True, err=True)
        # 可以在这里添加更详细的错误追踪信息
        # import traceback
        # click.echo(traceback.format_exc(), err=True)
        sys.exit(1)
    finally:
        # 无论成功失败，都恢复原始的sys.argv
        sys.argv = original_argv

if __name__ == '__main__':
    genomeasm()