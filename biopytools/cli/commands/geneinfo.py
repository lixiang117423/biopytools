"""
GFF3文件解析命令 | GFF3 File Parsing Command
"""

import click
import sys
from ...gff_utils.main import main as gff_main


@click.command(short_help = '从GFF文件中提取基因和转录本的信息',
               context_settings=dict(help_option_names=['-h', '--help'],max_content_width=120))
@click.option('--gff3', '-g',
              required=True,
              type=click.Path(exists=True),
              help='输入的GFF3文件路径 | Input GFF3 file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出的TSV文件路径 | Output TSV file path')
@click.option('--gene-type',
              default='gene',
              help='基因特征类型 (默认: gene) | Gene feature type (default: gene)')
@click.option('--transcript-types',
              multiple=True,
              default=['mRNA', 'transcript'],
              help='转录本特征类型 (默认: mRNA, transcript) | Transcript feature types (default: mRNA, transcript)')
def geneinfo(gff3, output, gene_type, transcript_types):
    """
    GFF3文件基因转录本信息提取工具
    
    从GFF3注释文件中提取基因和转录本的详细信息，
    输出为表格格式，便于后续分析使用。
    
    示例 | Examples:
    
    \b
    # 基本提取
    biopytools gff-utils -g genome.gff3 -o gene_info.tsv
    
    \b
    # 自定义特征类型
    biopytools gff-utils -g genome.gff3 -o gene_info.tsv \\
        --gene-type gene --transcript-types mRNA lnc_RNA
    
    \b
    # 提取所有转录本类型
    biopytools gff-utils -g genome.gff3 -o gene_info.tsv \\
        --transcript-types mRNA transcript lnc_RNA miRNA
    """
    
    # 构建参数列表传递给原始main函数
    args = ['geneinfo.py']
    
    # 必需参数
    args.extend(['-g', gff3])
    args.extend(['-o', output])
    
    # 可选参数（只在非默认值时添加）
    if gene_type != 'gene':
        args.extend(['--gene-type', gene_type])
    
    # transcript-types参数处理
    if transcript_types and transcript_types != ('mRNA', 'transcript'):
        args.extend(['--transcript-types'] + list(transcript_types))
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        gff_main()
    except SystemExit as e:
        # 处理程序正常退出
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv