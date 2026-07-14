"""
序列提取命令|Sequence extraction command
"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...seq_extract.main import main as seq_extract_main
        return seq_extract_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否帮助请求|Check if help request"""
    help_flags = {"-h", "--help"}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(非帮助模式)|Validate file exists (non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help="序列提取(seqkit封装,自动识别ID/ID文件/BED)|Sequence extraction (seqkit wrapper, auto-detect)",
    context_settings=dict(help_option_names=["-h", "--help"], max_content_width=120),
)
@click.option("-i", "--input",
              required=True,
              help="查询:单个ID、ID文件(一列)或BED文件(>=2列)|Query: single ID, ID file (1 column), or BED file (>=2 columns)")
@click.option("-s", "--sequence",
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value),
              help="目标序列FASTA文件|Target sequence FASTA file")
@click.option("-o", "--output",
              default=None,
              help="输出文件(默认自动推导:{query}.{subject}.fa)|Output file (default: auto-derived)")
@click.option("--bed",
              is_flag=True,
              help="强制BED模式(跳过自动检测)|Force BED mode (skip auto-detection)")
def seq_extract(input, sequence, output, bed):
    """序列提取工具(seqkit封装,自动识别ID/ID文件/BED)|Sequence extraction tool (seqkit wrapper, auto-detect ID/ID file/BED)

    示例|Examples: biopytools seq-extract -i gene.id.txt -s gene.fa -o gene.genomic.fa
    """
    seq_extract_main = _lazy_import_main()

    args = ["seq_extract.py"]
    args.extend(["-i", input])
    args.extend(["-s", sequence])
    if output:
        args.extend(["-o", output])
    if bed:
        args.append("--bed")

    original_argv = sys.argv
    sys.argv = args
    try:
        seq_extract_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
