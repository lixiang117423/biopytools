"""
🧬 转录组预测分析命令 | Transcriptome Prediction Analysis Command
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...transcriptome_prediction.main import main as transcriptome_prediction_main

# --- Placeholder for the original main function to make this snippet runnable ---
# --- In your project, delete this section and use the import statement above ---
# START: Placeholder for demonstration
def get_original_main_for_demo():
    def main_placeholder():
        print("--- 🚀 Original main function called (simulated) ---")
        print(f"Received sys.argv: {sys.argv}")
        # The original main() would create a parser and run the analysis.
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('-g', '--genome', required=True)
        input_group = parser.add_mutually_exclusive_group(required=True)
        input_group.add_argument('-r', '--rna-seq', nargs='+')
        input_group.add_argument('--samples-file')
        parser.add_argument('-o', '--output', required=True)
        parser.add_argument('-t', '--threads', type=int, default=88)
        # ... and all other args
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Analysis finished (simulated) ---")
    return main_placeholder
transcriptome_prediction_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']), short_help = "基于转录组的基因组转录本预测")
# --- Core Arguments ---
@click.option('--genome', '-g',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 基因组FASTA文件 | Genome FASTA file.')
@click.option('--rna-seq', '-r',
              multiple=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='📊 RNA-seq FASTQ文件 (可多个) | RNA-seq FASTQ files (multiple allowed).')
@click.option('--samples-file',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='📋 Trinity样本文件格式 | Trinity samples file format.')
@click.option('--output', '-o', 'output_dir',
              required=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='📂 输出目录 | Output directory.')
@click.option('--threads', '-t',
              type=int, default=88, show_default=True,
              help='🚀 线程数 | Number of threads.')

# --- Tool-specific Arguments (organized by comment blocks) ---
# --- HISAT2 Parameters ---
@click.option('--hisat2-min-intron', type=int, default=20, show_default=True, help='🎯 HISAT2最小内含子长度 | HISAT2 min intron length.')
@click.option('--hisat2-max-intron', type=int, default=500000, show_default=True, help='🎯 HISAT2最大内含子长度 | HISAT2 max intron length.')
@click.option('--dta/--no-dta', 'hisat2_dta', default=True, show_default=True, help='🧬 HISAT2 --dta 选项 (推荐用于StringTie) | HISAT2 --dta option.')

# --- StringTie Parameters ---
@click.option('--stringtie-min-length', type=int, default=200, show_default=True, help='🧩 StringTie最小转录本长度 | StringTie min transcript length.')

# --- Trinity Parameters ---
@click.option('--trinity-min-contig-length', type=int, default=200, show_default=True, help='🔗 Trinity最小contig长度 | Trinity min contig length.')
@click.option('--trinity-max-memory', default='20G', show_default=True, help='🧠 Trinity最大内存使用 | Trinity max memory usage.')
@click.option('--trinity-ss-lib-type', type=click.Choice(['FR', 'RF', 'F', 'R'], case_sensitive=False), help='⛓️ Trinity链特异性类型 | Trinity strand-specific library type.')

# --- TransDecoder Parameters ---
@click.option('--transdecoder-min-protein-len', type=int, default=100, show_default=True, help=' protein长度 | TransDecoder min protein length.')
def transcriptome_prediction(**kwargs):
    """🧬 转录组从头预测注释流程。
    
    一个集成的流程，使用HISAT2, StringTie, Trinity, PASA和TransDecoder
    从RNA-Seq数据中进行转录本的组装、注释和编码区预测。
    
    🌟 示例 | Examples:
    
    \b
    # 🎯 基本用法 (提供配对FASTQ文件)
    biopytools transcriptome-prediction -g g.fa -r s1_R1.fq -r s1_R2.fq -o results
    
    \b
    # 📋 使用样本文件作为输入
    biopytools transcriptome-prediction -g g.fa --samples-file samples.txt -o results
        
    \b
    # 🚀 高性能运行并调整参数
    biopytools transcriptome-prediction -g g.fa -r s1_R1.fq -r s1_R2.fq -o out -t 96 \\
        --trinity-max-memory 100G --stringtie-min-length 300
    """
    
    # 验证互斥参数 | Validate mutually exclusive arguments
    if not kwargs.get('rna_seq') and not kwargs.get('samples_file'):
        raise click.UsageError("❌ 必须提供 --rna-seq 或 --samples-file 中的一个 | Either --rna-seq or --samples-file must be provided.")
    if kwargs.get('rna_seq') and kwargs.get('samples_file'):
        raise click.UsageError("❌ --rna-seq 和 --samples-file 不能同时使用 | --rna-seq and --samples-file are mutually exclusive.")
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'transcriptome-prediction']
    
    # 手动处理所有参数以确保正确映射
    args.extend(['-g', kwargs['genome']])
    args.extend(['-o', kwargs['output_dir']])
    if kwargs.get('rna_seq'):
        args.append('-r')
        args.extend(kwargs['rna_seq'])
    if kwargs.get('samples_file'):
        args.extend(['--samples-file', kwargs['samples_file']])
    if kwargs.get('threads') != 88:
        args.extend(['-t', str(kwargs['threads'])])
    
    # HISAT2
    if kwargs.get('hisat2_min_intron') != 20: args.extend(['--hisat2-min-intron', str(kwargs['hisat2_min_intron'])])
    if kwargs.get('hisat2_max_intron') != 500000: args.extend(['--hisat2-max-intron', str(kwargs['hisat2_max_intron'])])
    if not kwargs.get('hisat2_dta'): args.append('--no-dta') # Handle inverted flag
    
    # StringTie
    if kwargs.get('stringtie_min_length') != 200: args.extend(['--stringtie-min-length', str(kwargs['stringtie_min_length'])])
    
    # Trinity
    if kwargs.get('trinity_min_contig_length') != 200: args.extend(['--trinity-min-contig-length', str(kwargs['trinity_min_contig_length'])])
    if kwargs.get('trinity_max_memory') != '20G': args.extend(['--trinity-max-memory', kwargs['trinity_max_memory']])
    if kwargs.get('trinity_ss_lib_type'): args.extend(['--trinity-ss-lib-type', kwargs['trinity_ss_lib_type']])
    
    # TransDecoder
    if kwargs.get('transdecoder_min_protein_len') != 100: args.extend(['--transdecoder-min-protein-len', str(kwargs['transdecoder_min_protein_len'])])
    
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        transcriptome_prediction_main()
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
    transcriptome_prediction()