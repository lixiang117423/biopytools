"""
🚀 GTX WGS分析命令 | GTX WGS Analysis Command 🚀
"""

import click
import sys
# In your actual project, you would use a relative import like this:
# from ...gtx_wgs.main import main as gtx_main

# --- Placeholder for demonstration purposes to make this script self-contained ---
# --- In your project, delete this section and use the import statement above ---
# START: Placeholder for demonstration
def get_original_main_for_demo():
    def main_placeholder():
        print("--- 🚀 Original main function called (simulated) ---")
        print(f"Received sys.argv: {sys.argv}")
        # Here, the original argparse would parse sys.argv and run the GTXAnalyzer.
        # We can simulate the parsing to verify the arguments.
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input-dir', required=True)
        parser.add_argument('-o', '--output-dir', required=True)
        parser.add_argument('-r', '--reference', required=True)
        parser.add_argument('-t', '--threads', type=int, default=88)
        parser.add_argument('--gtx-path', default='/share/apps/gtx/GTX.CAT_2.2.1/bin/gtx')
        parser.add_argument('--tmp-dir')
        parser.add_argument('--min-confidence', type=int, default=30)
        parser.add_argument('--min-base-quality', type=int, default=20)
        parser.add_argument('--ploidy', type=int, default=2)
        parser.add_argument('--pcr-indel-model', default='CONSERVATIVE')
        parser.add_argument('--read1-pattern', default='*_1.fq.gz')
        parser.add_argument('--read2-pattern', default='*_2.fq.gz')
        try:
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Analysis finished (simulated) ---")
    return main_placeholder
gtx_main = get_original_main_for_demo()
# END: Placeholder

@click.command(context_settings=dict(help_option_names=['-h', '--help']), short_help = "GTX WGS批处理分析流程")
@click.option('--input-dir', '-i',
              required=True,
              type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
              help='📂 输入目录 (包含clean FASTQ文件) | Input directory (containing clean FASTQ files).')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='📤 输出目录 | Output directory.')
@click.option('--reference', '-r',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='🧬 参考基因组文件 | Reference genome file.')
@click.option('--threads', '-t',
              type=int,
              default=88,
              show_default=True,
              help='🧵 线程数 | Number of threads.')
@click.option('--gtx-path',
              default='/share/apps/gtx/GTX.CAT_2.2.1/bin/gtx',
              show_default=True,
              help='💻 GTX程序路径 | GTX program path.')
@click.option('--tmp-dir',
              type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
              help='🗑️ 临时目录路径 (默认: output/tmp) | Temporary directory path.')
@click.option('--min-confidence',
              type=int,
              default=30,
              show_default=True,
              help='🎯 最小置信度阈值 | Minimum confidence threshold.')
@click.option('--min-base-quality',
              type=int,
              default=20,
              show_default=True,
              help='✨ 最小碱基质量阈值 | Minimum base quality threshold.')
@click.option('--ploidy',
              type=int,
              default=2,
              show_default=True,
              help='🧬 倍性 | Ploidy.')
@click.option('--pcr-indel-model',
              default='CONSERVATIVE',
              show_default=True,
              help='🔬 PCR indel模型 | PCR indel model.')
@click.option('--read1-pattern',
              default='*_1.fq.gz',
              show_default=True,
              help='📄 R1文件匹配模式 | R1 file pattern.')
@click.option('--read2-pattern',
              default='*_2.fq.gz',
              show_default=True,
              help='📄 R2文件匹配模式 | R2 file pattern.')
def gtx(input_dir, output_dir, reference, threads, gtx_path, tmp_dir,
        min_confidence, min_base_quality, ploidy, pcr_indel_model,
        read1_pattern, read2_pattern):
    """
    GTX WGS批处理分析流程.

    自动化执行从FASTQ到VCF的全基因组重测序分析，包括比对、
    排序、去重和变异检测。
    
    💡 示例 | Examples:
    
    \b
    # 🎯 基本运行 (必需参数)
    biopytools gtx -i ./clean_data -o ./gtx_results -r ref.fa
    
    \b
    # ⚡️ 使用更多线程和自定义GTX路径
    biopytools gtx -i ./data -o ./results -r hg38.fa -t 96 \\
        --gtx-path /opt/gtx/bin/gtx
    
    \b
    # 🔬 调整质量控制参数
    biopytools gtx -i ./data -o ./results -r ref.fa \\
        --min-confidence 25 --min-base-quality 15
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'gtx']
    
    # 必需参数 📌 | Required parameters
    args.extend(['-i', input_dir])
    args.extend(['-o', output_dir])
    args.extend(['-r', reference])
    
    # 可选参数 ✨ | Optional parameters
    if threads != 88:
        args.extend(['-t', str(threads)])
        
    if gtx_path != '/share/apps/gtx/GTX.CAT_2.2.1/bin/gtx':
        args.extend(['--gtx-path', gtx_path])
        
    if tmp_dir: # No default, so only add if provided
        args.extend(['--tmp-dir', tmp_dir])
        
    if min_confidence != 30:
        args.extend(['--min-confidence', str(min_confidence)])
        
    if min_base_quality != 20:
        args.extend(['--min-base-quality', str(min_base_quality)])
        
    if ploidy != 2:
        args.extend(['--ploidy', str(ploidy)])
        
    if pcr_indel_model != 'CONSERVATIVE':
        args.extend(['--pcr-indel-model', pcr_indel_model])
        
    if read1_pattern != '*_1.fq.gz':
        args.extend(['--read1-pattern', read1_pattern])
        
    if read2_pattern != '*_2.fq.gz':
        args.extend(['--read2-pattern', read2_pattern])

    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        gtx_main()
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
    # 模拟命令行调用，例如: python your_script.py -i . -o ./out -r ref.fa
    gtx()