"""
K-mer丰度分析命令 | K-mer Abundance Analysis Command
"""

import click
import sys

# 尝试导入，处理可能的依赖问题
try:
    from ...kmer_count.main import main as kmer_count_main
    IMPORT_SUCCESS = True
except ImportError as e:
    IMPORT_SUCCESS = False
    IMPORT_ERROR = str(e)
    def kmer_count_main():
        pass


@click.command(short_help='计算K-mer的丰度和滑窗分析',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='📁 输入文件目录 | Input files directory')
@click.option('--pattern', '-p',
              required=True,
              type=str,
              help='📁 文件模式，支持FASTQ和FASTA格式，如*_1.fq.gz、*.fa | File pattern, support FASTQ and FASTA formats')
@click.option('--kmer-lib', '-k',
              required=True,
              type=click.Path(exists=True),
              help='🧬 K-mer库文件(FASTA格式) | K-mer library file (FASTA format)')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='📂 输出目录 | Output directory')
@click.option('--bed-file', '-b',
              type=click.Path(exists=True),
              help='📋 BED文件路径 | BED file path')
@click.option('--kmer-size', '-m',
              default=51,
              type=int,
              help='📏 K-mer长度 | K-mer size (default: 51)')
@click.option('--hash-size', '-s',
              default='1000M',
              type=str,
              help='🗂️ 哈希表大小 | Hash table size (default: 1000M)')
@click.option('--threads', '-t',
              default=8,
              type=int,
              help='🧵 线程数 | Number of threads (default: 8)')
@click.option('--window-size', '-w',
              default=500000,
              type=int,
              help='🪟 滑动窗口大小bp | Sliding window size in bp (default: 500000)')
@click.option('--step-size',
              type=int,
              help='👣 滑动窗口步长bp (默认: window-size/5) | Sliding window step size in bp (default: window-size/5)')
@click.option('--canonical', '-C',
              is_flag=True,
              help='🔄 统计正向和反向互补链 | Count both forward and reverse complement')
@click.option('--keep-temp',
              is_flag=True,
              help='💾 保留临时文件 | Keep temporary files')
@click.option('--keep-binary',
              is_flag=True,
              help='🔢 保留0/1存在缺失矩阵 | Keep 0/1 presence/absence matrix')
@click.option('--jellyfish-path',
              default='jellyfish',
              type=str,
              help='🐙 Jellyfish程序路径 | Jellyfish program path (default: jellyfish)')
@click.option('--verbose', '-v',
              is_flag=True,
              help='📝 详细输出 | Verbose output')
def kmer_count(input, pattern, kmer_lib, output, bed_file, kmer_size, hash_size,
               threads, window_size, step_size, canonical, keep_temp, keep_binary,
               jellyfish_path, verbose):
    """
    🧬 K-mer丰度分析工具
    
    基于Jellyfish进行高效的k-mer计数分析，支持FASTQ和FASTA格式输入，
    可结合BED注释文件进行基因组区域特异性分析，支持滑动窗口分析。
    
    功能特点 | Features:
    - 高效的k-mer计数和丰度统计
    - 支持多种输入格式（FASTQ/FASTA）
    - BED文件注释支持
    - 滑动窗口分析
    - 多样本批处理
    - 0/1存在缺失矩阵生成
    
    输出文件 | Output Files:
    - kmer_abundance_matrix.tsv: 主要的丰度矩阵
    - sliding_window_analysis.tsv: 滑动窗口分析结果
    - kmer_binary_matrix.tsv: 0/1存在缺失矩阵
    - analysis_summary.txt: 分析统计报告
    - each_sample/: 每个样本的详细结果
    
    示例 | Examples:
    
    \b
    # 基本FASTQ分析
    biopytools kmer-count -i /data/fastq -p "*_1.fq.gz" \\
        -k kmers.fasta -o results/
    
    \b
    # 带BED注释的分析
    biopytools kmer-count -i ./samples -p "*_R1.fastq" \\
        -k kmers.fasta -b kmers.bed -o results/
    
    \b
    # FASTA格式输入
    biopytools kmer-count -i /data/fasta -p "*.fa" \\
        -k kmers.fasta -o results/
    
    \b
    # 滑动窗口分析
    biopytools kmer-count -i ./fastq -p "*_1.fq.gz" \\
        -k kmers.fasta -b regions.bed \\
        -w 500000 --step-size 100000 -o results/
    
    \b
    # 高线程数分析
    biopytools kmer-count -i /data/samples -p "*_R1.fastq.gz" \\
        -k target_kmers.fa -t 32 -s 4G -o analysis/
    
    \b
    # 正反链计数
    biopytools kmer-count -i ./reads -p "*_1.fq.gz" \\
        -k kmers.fasta -C --keep-binary -o results/
    
    \b
    # 自定义参数完整分析
    biopytools kmer-count -i /data/wgs -p "*_1.clean.fq.gz" \\
        -k target_sequences.fa -b annotation.bed \\
        -m 31 -s 2G -t 16 -w 1000000 --step-size 200000 \\
        --jellyfish-path /opt/jellyfish/bin/jellyfish \\
        --keep-temp --keep-binary --verbose -o comprehensive_analysis/
    
    文件格式说明 | File Format Notes:
    
    K-mer库文件 (FASTA):
    >kmer_id_1
    ATCGATCGATCG...
    >kmer_id_2  
    GCTAGCTAGCTA...
    
    BED文件格式:
    chr1    1000    2000    kmer_id_1
    chr1    3000    4000    kmer_id_2
    
    输出矩阵格式:
    chr    start    end    kmer    sample1    sample2    sample3
    chr1   1000     2000   ATCG... 15        23        8
    
    性能建议 | Performance Tips:
    - 大基因组推荐使用较大的哈希表大小 (-s 4G或更大)
    - 多样本分析时适当增加线程数
    - 使用SSD存储临时文件可显著提升速度
    - 对于大型数据集，考虑分批处理
    """
    
    # 构建参数列表传递给原始main函数
    args = ['kmer_count.py']
    
    # 必需参数
    args.extend(['-i', input])
    args.extend(['-p', pattern])
    args.extend(['-k', kmer_lib])
    args.extend(['-o', output])
    
    # 可选参数（只在非默认值时添加）
    if bed_file:
        args.extend(['-b', bed_file])
    
    if kmer_size != 51:
        args.extend(['-m', str(kmer_size)])
    
    if hash_size != '1000M':
        args.extend(['-s', hash_size])
    
    if threads != 8:
        args.extend(['-t', str(threads)])
    
    if window_size != 500000:
        args.extend(['-w', str(window_size)])
    
    if step_size:
        args.extend(['--step-size', str(step_size)])
    
    if canonical:
        args.append('-C')
    
    if keep_temp:
        args.append('--keep-temp')
    
    if keep_binary:
        args.append('--keep-binary')
    
    if jellyfish_path != 'jellyfish':
        args.extend(['--jellyfish-path', jellyfish_path])
    
    if verbose:
        args.append('-v')
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        kmer_count_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n分析被用户中断 | Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行错误 | Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv