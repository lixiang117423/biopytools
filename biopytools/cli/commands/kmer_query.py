"""
K-mer提取分析命令 | K-mer Extraction Analysis Command
"""

import click
import sys
from ...kmer_extractor.main import main as kmer_extract_main


@click.command(short_help='从FASTA/FASTQ文件中提取K-mer序列',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input-files', '-i',
              required=True,
              multiple=True,
              type=click.Path(exists=True),
              help='输入文件路径 (FASTA/FASTQ，支持压缩格式) | Input file paths (FASTA/FASTQ, supports compressed formats)')
@click.option('--output-dir', '-o',
              default='./kmer_output',
              type=click.Path(),
              help='输出目录 | Output directory (default: ./kmer_output)')
@click.option('--kmer-length', '-k',
              default=51,
              type=int,
              help='K-mer长度 (1-64) | K-mer length (1-64) (default: 51)')
@click.option('--threads', '-t',
              default=88,
              type=int,
              help='线程数 | Number of threads (default: 88)')
@click.option('--memory', '-m',
              default=880,
              type=int,
              help='内存限制(GB) | Memory limit (GB) (default: 880)')
@click.option('--file-type',
              type=click.Choice(['fasta', 'fastq']),
              help='文件类型 (自动检测如果未指定) | File type (auto-detect if not specified)')
@click.option('--fastq-pattern',
              type=str,
              help='FASTQ文件匹配模式 (例如: "*_1.clean.fq.gz") | FASTQ file matching pattern (e.g., "*_1.clean.fq.gz")')
@click.option('--no-canonical',
              is_flag=True,
              help='不使用canonical k-mer | Do not use canonical k-mers')
@click.option('--no-compress',
              is_flag=True,
              help='不压缩输出文件 | Do not compress output files')
@click.option('--output-bed',
              is_flag=True,
              help='输出BED格式文件 (仅适用于FASTA输入) | Output BED format file (only for FASTA input)')
@click.option('--no-keep-binary',
              is_flag=True,
              help='不保留二进制文件 (默认保留) | Do not keep binary files (keep by default)')
@click.option('--unikmer-path',
              default='unikmer',
              type=str,
              help='Unikmer软件路径 | Unikmer software path (default: unikmer)')
@click.option('--jellyfish-path',
              default='jellyfish',
              type=str,
              help='Jellyfish软件路径 | Jellyfish software path (default: jellyfish)')
@click.option('--jellyfish-hash-size',
              default='10000M',
              type=str,
              help='Jellyfish哈希表大小 | Jellyfish hash table size (default: 10000M)')
def kmer_query(input_files, output_dir, kmer_length, threads, memory, file_type,
                 fastq_pattern, no_canonical, no_compress, output_bed, no_keep_binary,
                 unikmer_path, jellyfish_path, jellyfish_hash_size):
    """
    K-mer提取工具 (FASTA使用unikmer，FASTQ使用jellyfish)
    
    从FASTA或FASTQ文件中高效提取K-mer序列，支持多种输入格式和输出选项。
    根据输入文件类型自动选择最优的提取工具：FASTA文件使用unikmer，FASTQ文件使用jellyfish。
    
    功能特点 | Features:
    - 支持FASTA和FASTQ格式输入（包括压缩文件）
    - 自动文件类型检测
    - FASTQ文件的配对模式匹配
    - 多线程并行处理
    - 灵活的输出格式选择
    - 内存使用优化
    
    输出文件 | Output Files:
    - <basename>.fasta: 提取的K-mer序列
    - <basename>.bed: K-mer位置注释文件（FASTA输入时可选）
    - <basename>.jf: Jellyfish二进制文件（FASTQ输入时可选保留）
    
    示例 | Examples:
    
    \b
    # 基本FASTA文件处理
    biopytools kmer-extract -i data.fasta -o results
    
    \b
    # 多个FASTQ文件处理
    biopytools kmer-extract -i sample1.fastq sample2.fastq -o results -k 31 -t 16
    
    \b
    # FASTQ配对文件处理
    biopytools kmer-extract -i /data/*.fastq.gz -o results \\
        --fastq-pattern "*_1.clean.fq.gz" -k 21
    
    \b
    # FASTA文件生成BED注释
    biopytools kmer-extract -i data.fasta -o results --output-bed -k 25
    
    \b
    # 多个FASTA文件批量处理
    biopytools kmer-extract --input-files file1.fa file2.fa \\
        -o results --threads 8 --memory 100 --output-bed
    
    \b
    # 高性能大数据处理
    biopytools kmer-extract -i large_genome.fa -o results \\
        -k 31 -t 128 -m 500 --unikmer-path /opt/unikmer/bin/unikmer
    
    \b
    # FASTQ文件自定义设置
    biopytools kmer-extract -i reads.fq.gz -o output \\
        --jellyfish-hash-size 20000M --no-canonical \\
        --jellyfish-path /usr/local/bin/jellyfish
    
    \b
    # 目录批量处理
    biopytools kmer-extract -i /data/genomes/ -o results \\
        --file-type fasta -k 51 --output-bed --no-compress
    
    文件格式说明 | File Format Notes:
    
    输入格式支持:
    - FASTA: .fasta, .fa, .fas (.gz压缩)
    - FASTQ: .fastq, .fq (.gz压缩)
    
    FASTQ配对模式示例:
    - "*_1.fq.gz" 匹配 sample_1.fq.gz 和 sample_2.fq.gz
    - "*_R1.fastq" 匹配 sample_R1.fastq 和 sample_R2.fastq
    
    输出FASTA格式:
    >kmer_1
    ATCGATCGATCGATCG...
    >kmer_2
    GCTAGCTAGCTAGCTA...
    
    BED文件格式 (FASTA输入):
    chr1    100    151    kmer_1
    chr1    200    251    kmer_2
    
    性能建议 | Performance Tips:
    - 大基因组建议增加内存限制 (-m 参数)
    - 使用SSD存储可显著提升I/O性能
    - 对于FASTQ文件，适当调整哈希表大小
    - 多文件处理时建议使用更多线程
    """
    
    # 构建参数列表传递给原始main函数
    args = ['kmer_query.py']
    
    # 必需参数
    args.extend(['-i'] + list(input_files))
    
    # 可选参数（只在非默认值时添加）
    if output_dir != './kmer_output':
        args.extend(['-o', output_dir])
    
    if kmer_length != 51:
        args.extend(['-k', str(kmer_length)])
    
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    if memory != 880:
        args.extend(['-m', str(memory)])
    
    if file_type:
        args.extend(['--file-type', file_type])
    
    if fastq_pattern:
        args.extend(['--fastq-pattern', fastq_pattern])
    
    if no_canonical:
        args.append('--no-canonical')
    
    if no_compress:
        args.append('--no-compress')
    
    if output_bed:
        args.append('--output-bed')
    
    if no_keep_binary:
        args.append('--no-keep-binary')
    
    if unikmer_path != 'unikmer':
        args.extend(['--unikmer-path', unikmer_path])
    
    if jellyfish_path != 'jellyfish':
        args.extend(['--jellyfish-path', jellyfish_path])
    
    if jellyfish_hash_size != '10000M':
        args.extend(['--jellyfish-hash-size', jellyfish_hash_size])
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        kmer_extract_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n提取流程被用户中断 | Extraction pipeline interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行错误 | Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv