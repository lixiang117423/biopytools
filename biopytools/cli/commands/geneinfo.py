"""
GFF3文件解析命令 | GFF3 File Parsing Command
优化版本：使用懒加载解决响应速度问题
"""

import click
import sys
import os


def _lazy_import_gff_main():
    """懒加载gff main函数 | Lazy load gff main function"""
    try:
        from ...gff_utils.main import main as gff_main
        return gff_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_gff3_file(file_path):
    """验证GFF3文件是否存在（仅在非帮助模式下）| Validate GFF3 file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"GFF3文件不存在 | GFF3 file does not exist: {file_path}")
    return file_path


def _validate_output_path(file_path):
    """验证输出路径（仅在非帮助模式下）| Validate output path (only in non-help mode)"""
    if not _is_help_request():
        # 输出文件可以不存在，但其父目录应该存在或可创建
        parent_dir = os.path.dirname(os.path.abspath(file_path))
        if parent_dir and not os.path.exists(parent_dir):
            try:
                os.makedirs(parent_dir, exist_ok=True)
            except OSError as e:
                raise click.BadParameter(f"无法创建输出目录 | Cannot create output directory: {parent_dir}, Error: {e}")
    return file_path


@click.command(short_help='从GFF文件中提取基因和转录本的信息',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--gff3', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_gff3_file(value) if value else None,
              help='输入的GFF3文件路径 | Input GFF3 file path')
@click.option('--output', '-o',
              required=True,
              callback=lambda ctx, param, value: _validate_output_path(value) if value else None,
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
    biopytools geneinfo -g genome.gff3 -o gene_info.tsv
    
    \b
    # 自定义特征类型
    biopytools geneinfo -g genome.gff3 -o gene_info.tsv \\
        --gene-type gene --transcript-types mRNA lnc_RNA
    
    \b
    # 提取所有转录本类型
    biopytools geneinfo -g genome.gff3 -o gene_info.tsv \\
        --transcript-types mRNA transcript lnc_RNA miRNA
    """
    
    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    gff_main = _lazy_import_gff_main()
    
    # 构建参数列表传递给原始main函数
    args = ['geneinfo.py']
    
    # 必需参数
    args.extend(['-g', gff3])
    args.extend(['-o', output])
    
    # 可选参数（只在非默认值时添加）
    if gene_type != 'gene':
        args.extend(['--gene-type', gene_type])
    
    # transcript-types参数处理
    # 注意：比较tuple时要考虑顺序和内容
    default_transcript_types = ('mRNA', 'transcript')
    if transcript_types and tuple(transcript_types) != default_transcript_types:
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