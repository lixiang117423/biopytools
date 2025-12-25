"""
Genome Analysis Command
基因组分析命令
"""

import click
import sys
import os


def _lazy_import_genome_analysis_main():
    """懒加载genome_analysis main函数 | Lazy load genome_analysis main function"""
    try:
        from ...genome_analysis.main import main as genome_analysis_main
        return genome_analysis_main
    except ImportError as e:
        click.echo(f"[ERROR] 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_dir_exists(dir_path):
    """验证目录是否存在（仅在非帮助模式下）| Validate directory exists (only in non-help mode)"""
    if _is_help_request():
        return dir_path
    if not os.path.exists(dir_path):
        raise click.BadParameter(f"目录不存在 | Directory does not exist: {dir_path}")
    if not os.path.isdir(dir_path):
        raise click.BadParameter(f"路径不是目录 | Path is not a directory: {dir_path}")
    return dir_path


@click.command(short_help="GenomeScope2基因组评估工具",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input-dir', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_dir_exists(value) if value else None,
              help='[DIR] 输入FASTQ文件目录 | Input FASTQ directory')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(),
              help='[DIR] 输出目录 | Output directory')
@click.option('--read-length', '-l',
              type=int,
              default=150,
              help='[INT] 测序读长 (默认: 150) | Read length (default: 150)')
@click.option('--kmer-size', '-k',
              type=int,
              default=21,
              help='[INT] K-mer大小 (默认: 21) | K-mer size (default: 21)')
@click.option('--threads', '-t',
              type=int,
              default=64,
              help='[INT] 线程数 (默认: 64) | Number of threads (default: 64)')
@click.option('--hash-size', '-s',
              default='10G',
              help='[STR] Jellyfish哈希表大小 (默认: 10G) | Jellyfish hash size (default: 10G)')
@click.option('--max-kmer-cov', '-c',
              type=int,
              default=1000,
              help='[INT] 最大k-mer覆盖度 (默认: 1000) | Max k-mer coverage (default: 1000)')
@click.option('--genomescope-r',
              default='~/software/scripts/genomescope.R',
              help='[FILE] GenomeScope R脚本路径 | GenomeScope R script path')
@click.option('--run-smudgeplot',
              is_flag=True,
              default=False,
              help='[FLAG] 运行Smudgeplot倍性分析 | Run Smudgeplot ploidy analysis')
@click.option('--fastk-table',
              default='',
              help='[FILE] FastK表文件路径 (如已存在) | FastK table file path (if exists)')
@click.option('--fastk-memory',
              default='16G',
              help='[STR] FastK内存大小 (默认: 16G) | FastK memory size (default: 16G)')
def genome_analysis(input_dir, output_dir, read_length, kmer_size, threads,
                   hash_size, max_kmer_cov, genomescope_r, run_smudgeplot,
                   fastk_table, fastk_memory):
    """
    GenomeScope2基因组评估工具

    使用Jellyfish和GenomeScope2评估基因组大小、杂合度等特征，
    并自动提取k-mer coverage用于后续分析。

    可选运行Smudgeplot进行倍性分析。

    示例 | Examples:

    \b
    # 基本用法 - 仅GenomeScope分析
    biopytools genome-analysis \\
        -i fastq_dir \\
        -o output_dir \\
        -l 150

    \b
    # 指定k-mer大小和线程数
    biopytools genome-analysis \\
        -i fastq_dir \\
        -o output_dir \\
        -l 150 \\
        -k 21 \\
        -t 32

    \b
    # 运行完整的GenomeScope + Smudgeplot分析
    biopytools genome-analysis \\
        -i fastq_dir \\
        -o output_dir \\
        -l 150 \\
        --run-smudgeplot

    \b
    # 使用已存在的FastK表运行Smudgeplot
    biopytools genome-analysis \\
        -i fastq_dir \\
        -o output_dir \\
        -l 150 \\
        --run-smudgeplot \\
        --fastk-table /path/to/fastk_table

    输出说明 | Output Description:

    工具会生成以下文件：

    \b
    - genome_analysis.jf          - Jellyfish计数结果
    - genome_analysis.histo       - K-mer直方图
    - genome_analysis_genomescope_output/  - GenomeScope结果目录
      - summary.txt              - 基因组特征汇总
      - model.txt                - 模型参数（包含kcov值）
      - various plots             - 可视化图表

    Smudgeplot输出 (使用--run-smudgeplot时):

    \b
    - fastk_table                 - FastK k-mer数据库
    - smudgeplot_result_kmerpairs.smu  - K-mer对文件
    - smudgeplot_result_*.pdf     - 倍性分析图表

    提取的参数 | Extracted Parameters:

    \b
    - kcov (k-mer coverage)       - 从GenomeScope model.txt提取
    - 基因组大小                  - 估算的基因组大小
    - 杂合度                      - 杂合度水平
    - 重复序列比例                - 重复序列占比
    - 倍性推断                    - Smudgeplot推断的倍性

    注意事项 | Notes:

    \b
    - K-mer大小建议使用21 (默认值)
    - 哈希表大小应根据基因组大小调整
    - Smudgeplot需要FastK工具和额外的计算时间
    - 输入目录中的所有.fastq/.fq/.fastq.gz/.fq.gz文件都会被使用
    """

    # 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    genome_analysis_main = _lazy_import_genome_analysis_main()

    # 构建参数列表传递给原始main函数 | Build argument list for original main function
    args = ['genome_analysis.py']

    # 必需参数 | Required parameters
    args.extend(['-i', input_dir])
    args.extend(['-o', output_dir])
    args.extend(['-l', str(read_length)])

    # 可选参数（只在非默认值时添加，减少命令行长度）| Optional parameters (add only when non-default)
    if kmer_size != 21:
        args.extend(['-k', str(kmer_size)])

    if threads != 64:
        args.extend(['-t', str(threads)])

    if hash_size != '10G':
        args.extend(['-s', hash_size])

    if max_kmer_cov != 1000:
        args.extend(['-c', str(max_kmer_cov)])

    if genomescope_r != '~/software/scripts/genomescope.R':
        args.extend(['--genomescope-r', genomescope_r])

    if run_smudgeplot:
        args.append('--run-smudgeplot')

    if fastk_table:
        args.extend(['--fastk-table', fastk_table])

    if fastk_memory != '16G':
        args.extend(['--fastk-memory', fastk_memory])

    # 保存并恢复sys.argv | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始的main函数 | Call original main function
        genome_analysis_main()
    except SystemExit as e:
        # 处理程序正常退出 | Handle normal program exit
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n[WARNING] 用户中断操作 | User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"[ERROR] 运行错误 | Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
