"""
基因序列提取命令 | Gene Sequence Extraction Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_gene_main():
    """懒加载基因提取main函数 | Lazy load gene extraction main function"""
    try:
        from ...parse_gene_seq.main import main as gene_main
        return gene_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"❌ 文件不存在 | File does not exist: {file_path}")
    return file_path


@click.command(short_help='🧬 基因序列提取工具：从基因组和GFF文件中提取基因序列',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='🧬 基因组FASTA文件路径 | Genome FASTA file path')
@click.option('--gff', '-f',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📄 GFF注释文件路径 | GFF annotation file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='💾 输出FASTA文件路径 | Output FASTA file path')
@click.option('--feature-type',
              default='gene',
              help='🎯 要提取的特征类型 | Feature type to extract (default: gene)')
@click.option('--min-length',
              type=int,
              default=0,
              help='📏 最小基因长度过滤 | Minimum gene length filter (default: 0)')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('--line-width',
              type=int,
              default=60,
              help='📐 FASTA序列行宽度 | FASTA sequence line width (default: 60)')
@click.option('--verbose', '-v',
              is_flag=True,
              help='🔍 显示详细信息 | Show verbose output')
def parse_gene_dna(genome, gff, output, feature_type, min_length, threads, line_width, verbose):
    """
    🧬 基因序列提取工具 | Gene Sequence Extraction Tool
    
    从基因组FASTA文件和GFF注释文件中提取指定类型的基因序列，
    支持多种特征类型提取和灵活的过滤条件。
    
    💡 示例 | Examples:
    
    \b
    # 🎯 提取所有基因序列
    biopytools parse-gene-dna -g genome.fasta -f annotation.gff -o genes.fasta
    
    \b
    # 🔬 提取CDS序列
    biopytools parse-gene-dna -g genome.fa -f genes.gff3 -o cds.fa --feature-type CDS
    
    \b
    # 📏 过滤短序列
    biopytools parse-gene-dna -g genome.fasta -f annotation.gff -o genes.fasta \\
        --min-length 300 -v
    """
    
    # 🚀 懒加载
    gene_main = _lazy_import_gene_main()
    
    # 构建参数列表
    args = ['parse_gene_dna.py']
    args.extend(['-g', genome])
    args.extend(['-f', gff])
    args.extend(['-o', output])
    
    if feature_type != 'gene':
        args.extend(['--feature-type', feature_type])
    if min_length != 0:
        args.extend(['--min-length', str(min_length)])
    if threads != 88:
        args.extend(['-t', str(threads)])
    if line_width != 60:
        args.extend(['--line-width', str(line_width)])
    if verbose:
        args.append('-v')
    
    original_argv = sys.argv
    sys.argv = args
    
    try:
        gene_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 基因序列提取被用户中断", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 基因序列提取失败: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv