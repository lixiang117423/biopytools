"""
VCF序列提取命令|VCF Sequence Extraction Command
"""

import click
import sys


def _lazy_import_vcf_sequence_main():
    """延迟加载vcf_sequence_toolkit主函数|Lazy load vcf_sequence_toolkit main function"""
    try:
        from ...vcf_sequence.main import main as vcf_sequence_main
        return vcf_sequence_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


@click.command(
    short_help='VCF序列提取工具|VCF sequence extraction tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--vcf', '-v',
              required=True,
              type=click.Path(exists=True),
              help='输入VCF文件路径|Input VCF file path')
@click.option('--genome', '-g',
              required=True,
              type=click.Path(exists=True),
              help='基因组FASTA文件路径|Genome FASTA file path')
@click.option('--chrom', '-c',
              required=True,
              type=str,
              help='染色体名称|Chromosome name')
@click.option('--start', '-s',
              required=True,
              type=int,
              help='起始位置|Start position')
@click.option('--end', '-e',
              required=True,
              type=int,
              help='结束位置|End position')
@click.option('--output-dir', '-o',
              default='./sequence_output',
              type=click.Path(),
              show_default=True,
              help='输出目录|Output directory')
@click.option('--format',
              type=click.Choice(['tab', 'fasta', 'csv']),
              default='tab',
              show_default=True,
              help='输出格式|Output format')
@click.option('--second-allele',
              is_flag=True,
              help='使用第二个等位基因|Use second allele')
@click.option('--no-reference',
              is_flag=True,
              help='不包含参考序列|Do not include reference sequence')
@click.option('--min-qual',
              type=int,
              help='最小质量过滤|Minimum quality filter')
@click.option('--samples',
              type=str,
              help='指定样本|Specify samples')
@click.option('--exclude-samples',
              type=str,
              help='排除样本|Exclude samples')
def vcf_sequence(vcf, genome, chrom, start, end, output_dir, format, second_allele,
                 no_reference, min_qual, samples, exclude_samples):
    """
    VCF序列提取工具|VCF Sequence Extraction Tool

    从VCF和基因组文件中提取特定区间的序列|Extract sequences from VCF and genome files for specific regions

    示例|Examples: biopytools vcf-sequence -v variants.vcf -g genome.fa -c chr1 -s 1000000 -e 1001000 -o seq_output
    """

    # 延迟加载|Lazy loading
    vcf_sequence_main = _lazy_import_vcf_sequence_main()

    # 构建参数列表|Build argument list
    args = ['vcf_sequence.py']

    # 必需参数|Required parameters
    args.extend(['-v', vcf])
    args.extend(['-g', genome])
    args.extend(['-c', chrom])
    args.extend(['-s', str(start)])
    args.extend(['-e', str(end)])

    # 可选参数（只在非默认值时添加）|Optional parameters (add only when non-default)
    if output_dir != './sequence_output':
        args.extend(['-o', output_dir])

    if format != 'tab':
        args.extend(['--format', format])

    if min_qual is not None:
        args.extend(['--min-qual', str(min_qual)])

    if samples:
        args.extend(['--samples', samples])

    if exclude_samples:
        args.extend(['--exclude-samples', exclude_samples])

    # 标志参数|Flag parameters
    if second_allele:
        args.append('--second-allele')

    if no_reference:
        args.append('--no-reference')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        result = vcf_sequence_main()
        sys.exit(result if result is not None else 0)
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(130)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
