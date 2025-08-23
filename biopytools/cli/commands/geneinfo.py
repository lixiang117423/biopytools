"""
🧬 GFF3工具命令 | GFF3 Tools Command
"""

import click
import sys
# In your actual project, you would import the original main function:
# from ...gff_tools.main import main as gff_main

# --- Placeholder for the original main function to make this snippet runnable ---
# --- In your project, delete this section and use the import statement above ---
# START: Placeholder for demonstration
def get_original_main_for_demo():
    def main_placeholder():
        print("--- 🚀 Original main function called (simulated) ---")
        print(f"Received sys.argv: {sys.argv}")
        # In a real run, argparse would parse sys.argv and the GFFAnalyzer would execute.
        # Example of how it would be parsed:
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('--gff3', '-g', required=True)
        parser.add_argument('--output', '-o', required=True)
        parser.add_argument('--gene-type', default='gene')
        parser.add_argument('--transcript-types', nargs='+', default=['mRNA', 'transcript'])
        try:
            # We don't run the full logic, just parse to show it works
            args = parser.parse_args()
            print(f"Argparse would have parsed arguments as: {args}")
        except Exception as e:
            print(f"Argparse simulation failed: {e}")
        
        print("--- ✅ Extraction finished (simulated) ---")
    return main_placeholder
gff_main = get_original_main_for_demo()
# END: Placeholder

@click.command(short_help = "从GFF3文件提取整合的基因和转录本信息")
@click.option('--gff3', '-g',
              required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help='📄 输入的GFF3文件路径 | Input GFF3 file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(dir_okay=False, resolve_path=True),
              help='💾 输出的TSV文件路径 | Output TSV file path')
@click.option('--gene-type',
              default='gene',
              show_default=True,
              help='🏷️ 基因特征类型 | Gene feature type')
@click.option('--transcript-types',
              multiple=True,
              default=['mRNA', 'transcript'],
              show_default=True,
              help='📜 转录本特征类型 (可多次使用) | Transcript feature types (use multiple times)')
def geneinfo(gff3, output, gene_type, transcript_types):
    """
    从GFF3文件提取整合的基因和转录本信息.

    该工具解析GFF3文件，将每个转录本与其对应的基因信息关联起来，
    并将结果输出到一个易于使用的TSV表格中。
    
    示例 | Examples:
    
    \b
    # 🎯 基本用法
    biopytools gff-tools -g input.gff3 -o gene_info.tsv
    
    \b
    # 📜 自定义转录本类型 (例如，处理lnc_RNA和tRNA)
    biopytools gff-tools -g annotation.gff3 -o custom_info.tsv \\
        --transcript-types mRNA \\
        --transcript-types lnc_RNA \\
        --transcript-types tRNA
    """
    
    # 构建参数列表以传递给原始的main函数 🔄 | Build argument list for original main function
    args = ['biopytools', 'gff-tools']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-g', gff3])
    args.extend(['-o', output])
    
    # 可选参数 ⚙️ | Optional parameters
    if gene_type != 'gene':
        args.extend(['--gene-type', gene_type])
        
    # 处理多值参数 (nargs='+') | Handle multiple value parameter (nargs='+')
    # click's 'multiple=True' collects values into a tuple. We need to add the flag
    # once, followed by all values for argparse's 'nargs='+' to parse correctly.
    # We compare sorted lists/tuples to handle order-insensitivity.
    if sorted(transcript_types) != sorted(['mRNA', 'transcript']):
        args.append('--transcript-types')
        args.extend(transcript_types)
        
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        gff_main()
    except SystemExit as e:
        # 处理程序正常退出 (如 argparse 的 --help) ✅ | Handle normal program exit
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
    # 模拟命令行调用，例如: python your_script.py -g file.gff3 -o out.tsv
    geneinfo()