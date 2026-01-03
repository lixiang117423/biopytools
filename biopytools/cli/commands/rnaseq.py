"""
RNA-seq分析命令|RNA-seq Analysis Command
"""

import click
import sys


@click.command(short_help='RNA-seq分析流程|RNA-seq analysis pipeline',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--genome', '-g',
              required=True,
              type=click.Path(exists=True),
              help='基因组FASTA文件路径|Genome FASTA file path')
@click.option('--gtf', '-f',
              required=True,
              type=click.Path(exists=True),
              help='基因注释GTF文件路径|Gene annotation GTF file path')
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='输入FASTQ文件目录或样本信息文件|Input FASTQ file directory or sample information file')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--pattern', '-p',
              type=str,
              help='FASTQ文件命名模式 (如: "*.R1.fastq.gz" 或 "*_1.fq.gz")|'
                   'Fastq file naming pattern (e.g., "*.R1.fastq.gz" or "*_1.fq.gz"), * represents sample name')
@click.option('--remove', '-r',
              default='no',
              type=click.Choice(['yes', 'y', 'no', 'n']),
              help='处理后删除BAM文件 (默认: no)|Remove BAM files after processing (default: no)')
@click.option('--verbose', '-v',
              count=True,
              help='增加输出详细程度 (-v, -vv, -vvv)|Increase output verbosity')
@click.option('--quiet',
              is_flag=True,
              help='静默模式，仅输出错误信息|Quiet mode, only output errors')
@click.option('--log-file',
              type=str,
              help='日志文件路径|Log file path')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
              default='INFO',
              help='日志级别|Log level')
@click.option('--dry-run',
              is_flag=True,
              help='试运行模式，不实际执行|Dry run mode, no actual execution')
@click.option('--threads', '-t',
              type=int,
              default=8,
              help='线程数 (默认: 8)|Number of threads (default: 8)')
@click.option('--sample-timeout',
              type=int,
              default=21600,
              help='单个样本处理超时时间（秒），默认6小时|Sample processing timeout in seconds (default: 21600 = 6 hours)')
def rnaseq(genome, gtf, input, output, pattern, remove, verbose, quiet,
           log_file, log_level, dry_run, threads, sample_timeout):
    """
    RNA-seq分析流程：HISAT2 + StringTie|RNA-seq Analysis Pipeline: HISAT2 + StringTie

    完整的RNA-seq数据分析流程，包括HISAT2比对和StringTie定量，计算FPKM和TPM值

    使用示例|Examples:

    \b
    # 基本用法
    biopytools rnaseq -g genome.fa -f genes.gtf \\
        -i /data/fastq/ -o rnaseq_results

    \b
    # 指定FASTQ文件命名模式
    biopytools rnaseq --genome reference.fasta --gtf annotation.gtf \\
        --input ./samples/ --output ./analysis/ \\
        --pattern "*.R1.fastq.gz"

    \b
    # 处理后删除BAM文件节省空间
    biopytools rnaseq -g genome.fa -f genes.gtf \\
        -i /data/rna_samples/ -o results/ \\
        -t 32 --remove yes

    \b
    # 设置样本处理超时时间（3小时）
    biopytools rnaseq -g genome.fa -f genes.gtf \\
        -i /data/rna_samples/ -o results/ \\
        --sample-timeout 10800
    """

    # Lazy loading: 只在调用时导入主模块|Import main module only when called
    from ...rnaseq.main import main as rnaseq_main

    # 构建参数列表|Build argument list
    args = ['rnaseq.py']

    # 必需参数|Required arguments
    args.extend(['-g', genome])
    args.extend(['-f', gtf])
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数|Optional arguments
    if pattern:
        args.extend(['-p', pattern])

    if remove != 'no':
        args.extend(['-r', remove])

    if threads != 8:
        args.extend(['-t', str(threads)])

    if sample_timeout != 21600:
        args.extend(['--sample-timeout', str(sample_timeout)])

    # 日志选项|Logging options
    if verbose > 0:
        args.extend(['-v'] * verbose)

    if quiet:
        args.append('--quiet')

    if log_file:
        args.extend(['--log-file', log_file])

    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    # 高级选项|Advanced options
    if dry_run:
        args.append('--dry-run')

    # 临时修改sys.argv|Temporarily modify sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        rnaseq_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断操作|Operation interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
