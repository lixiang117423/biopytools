"""
最长转录本提取CLI包装器|Longest mRNA Extraction CLI Wrapper
"""

import click
import sys
import os


def _lazy_import_longest_mrna_main():
    """延迟加载最长转录本提取主函数|Lazy load longest mRNA main function"""
    try:
        from ...longestmrna.main import main as longest_mrna_main
        return longest_mrna_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='最长转录本提取|Longest mRNA Extraction',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--genome', '-g',
              required=True,
              type=click.Path(exists=True),
              help='输入基因组FASTA文件|Input genome FASTA file')
@click.option('--gff3', '-f',
              required=True,
              type=click.Path(exists=True),
              help='输入GFF3注释文件|Input GFF3 annotation file')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出FASTA文件|Output FASTA file')
@click.option('--gene-info',
              type=click.Path(),
              help='基因信息输出文件|Gene info output file')
@click.option('--cds-output',
              type=click.Path(),
              help='CDS核苷酸序列输出文件|CDS nucleotide sequence output file')
def longestmrna(genome, gff3, output, gene_info, cds_output):
    """
    最长转录本提取流程|Longest mRNA Extraction Pipeline

    从GFF3注释文件和基因组FASTA文件中提取每个基因的最长mRNA转录本对应的CDS序列
    Extract longest mRNA CDS sequences for each gene from GFF3 annotation and genome FASTA

    示例|Examples: biopytools longest-mrna -g genome.fasta -f annotation.gff3 -o longest_proteins.fasta
    """

    # 延迟加载|Lazy load
    longest_mrna_main = _lazy_import_longest_mrna_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['longestmrna.py']

    # 必需参数|Required parameters
    args.extend(['-g', genome])
    args.extend(['-f', gff3])
    args.extend(['-o', output])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if gene_info:
        args.extend(['--gene-info', gene_info])
    if cds_output:
        args.extend(['--cds-output', cds_output])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        longest_mrna_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Extraction pipeline interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
