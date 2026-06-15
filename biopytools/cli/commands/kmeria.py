"""
K-mer GWAS分析命令|K-mer GWAS Analysis Command
"""

import click
import sys


@click.group(
    name='kmeria',
    short_help='K-mer GWAS全流程分析工具|K-mer GWAS Complete Pipeline Tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.pass_context
def kmeria(ctx):
    """
    K-mer GWAS全流程分析工具|K-mer GWAS Complete Pipeline Tool

    基于k-mer的全基因组关联分析工具，支持完整流程和各步骤独立运行|K-mer based GWAS tool supporting complete pipeline and individual steps
    """
    pass


@kmeria.command('pipeline')
@click.option('-i', '--fastq-dir', required=True,
              help='FASTQ文件目录|FASTQ files directory')
@click.option('--sample', 'samples', required=True,
              help='样本列表文件|Sample list file')
@click.option('-d', '--depth-file', required=True,
              help='测序深度文件|Sequencing depth file')
@click.option('-p', '--pheno-file', required=True,
              help='表型文件|Phenotype file')
@click.option('-o', '--output-dir', default='./kmeria_results',
              help='输出目录|Output directory')
@click.option('-f', '--force', is_flag=True, default=False,
              help='强制重新运行所有步骤|Force re-run all steps even if output exists')
@click.option('-k', '--kmer-size', default=31,
              help='K-mer大小|K-mer size (default: 31)')
@click.option('--max-abund', default=1000,
              help='最大丰度|Maximum abundance (default: 1000)')
@click.option('--missing-ratio', default=0.05,
              help='缺失率|Missing ratio (default: 0.05)')
@click.option('--ploidy', default=2,
              help='倍性|Ploidy (default: 2)')
@click.option('--step', type=click.Choice(['count', 'kctm', 'filter', 'm2b', 'asso']),
              help='从指定步骤开始|Start from specified step')
@click.option('-t', '--threads', default=12,
              help='线程数|Thread count (default: 12)')
@click.option('--batch-size', default=4,
              help='批处理大小|Batch size (default: 4)')
@click.option('--pheno-col', default=1,
              help='表型列|Phenotype column (default: 1)')
@click.option('--kinship-file',
              help='亲缘关系矩阵|Kinship matrix file')
@click.option('--covar-file',
              help='协变量文件|Covariate file')
@click.option('--enable-qc', is_flag=True, default=True,
              help='启用质控|Enable QC (default: True)')
@click.option('--enable-visualization', is_flag=True, default=True,
              help='启用可视化|Enable visualization (default: True)')
@click.option('--enable-annotation', is_flag=True,
              help='启用k-mer注释|Enable k-mer annotation')
@click.option('--genome-file',
              help='参考基因组|Reference genome (for annotation)')
@click.option('--gff-file',
              help='GFF注释文件|GFF annotation file (for annotation)')
@click.option('--sample-ratio', type=float, default=0.1,
              help='高p值位点抽样比例，用于减少绘图点数 (默认: 0.1 = 10%)|Sampling ratio for high p-value loci to reduce plot points (default: 0.1 = 10%)')
@click.option('--window-size', type=int, default=200000,
              help='基因查找窗口大小，单位bp (默认: 200000 = 200kb)|Gene search window size in bp (default: 200000 = 200kb)')
@click.option('--alignment-tool', type=click.Choice(['bwa', 'blast']), default='bwa',
              help='Post-GWAS比对工具选择 (默认: bwa)|Alignment tool for Post-GWAS (default: bwa)')
@click.option('--bwa-k', type=int, default=9,
              help='BWA mem -k 参数，最小种子长度 (默认: 9)|BWA mem -k parameter, minimum seed length (default: 9)')
@click.option('--bwa-t-min-score', type=int, default=10,
              help='BWA mem -T 参数，最小输出分数 (默认: 10)|BWA mem -T parameter, minimum score to output (default: 10)')
@click.option('--as-ratio', type=float, default=0.95,
              help='BWA AS过滤阈值，保留AS >= 最高AS * as_ratio的所有比对 (默认: 0.95)|BWA AS filtering threshold, keep alignments with AS >= max_AS * ratio (default: 0.95)')
@click.option('--log-file',
              help='日志文件|Log file')
def pipeline(fastq_dir, samples, depth_file, pheno_file, output_dir,
            kmer_size, max_abund, missing_ratio, ploidy,
            step, threads, batch_size, pheno_col, kinship_file, covar_file,
            enable_qc, enable_visualization, enable_annotation,
            genome_file, gff_file, sample_ratio, window_size,
            alignment_tool, bwa_k, bwa_t_min_score, as_ratio, log_file, force):
    """
    K-mer GWAS完整分析流程|K-mer GWAS Complete Pipeline

    示例|Example: biopytools kmeria pipeline -i fastq/ --samples samples.txt -d depth.txt -p pheno.txt -o results/
    """
    from biopytools.kmeria.main import main as kmeria_main

    args = ['kmeria.py', 'pipeline',
            '-i', fastq_dir,
            '--samples', samples,
            '-d', depth_file,
            '-p', pheno_file,
            '-o', output_dir,
            '-k', str(kmer_size),
            '--max-abund', str(max_abund),
            '--missing-ratio', str(missing_ratio),
            '--ploidy', str(ploidy),
            '-t', str(threads),
            '--batch-size', str(batch_size),
            '--pheno-col', str(pheno_col),
            '--sample-ratio', str(sample_ratio),
            '--window-size', str(window_size),
            '--alignment-tool', alignment_tool,
            '--bwa-k', str(bwa_k),
            '--bwa-t-min-score', str(bwa_t_min_score),
            '--as-ratio', str(as_ratio)]

    if step:
        args.extend(['--step', step])
    if force:
        args.append('--force')
    if kinship_file:
        args.extend(['--kinship-file', kinship_file])
    if covar_file:
        args.extend(['--covar-file', covar_file])

    # 可选功能参数（main.py中默认都是True，所以不需要显式传递）
    # Optional features (all default to True in main.py, so no need to pass explicitly)
    if enable_annotation:
        args.append('--enable-annotation')
    if genome_file:
        args.extend(['--genome-file', genome_file])
    if gff_file:
        args.extend(['--gff-file', gff_file])
    if log_file:
        args.extend(['--log-file', log_file])

    original_argv = sys.argv
    sys.argv = args

    try:
        kmeria_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv


@kmeria.command('count')
@click.option('-i', '--fastq-dir', required=True,
              help='FASTQ文件目录|FASTQ files directory')
@click.option('--sample', 'samples', required=True,
              help='样本列表文件|Sample list file')
@click.option('-o', '--output-dir', default='./01_kmer_counts',
              help='输出目录|Output directory')
@click.option('-k', '--kmer-size', default=31,
              help='K-mer大小|K-mer size (default: 31)')
@click.option('-t', '--threads', default=12,
              help='线程数|Thread count (default: 12)')
@click.option('-b', '--batch-size', default=4,
              help='批处理大小|Batch size (default: 4)')
@click.option('-C', '--count-separate-strands', is_flag=True,
              help='分别计数链|Count strands separately')
@click.option('-T', '--text-output', is_flag=True,
              help='文本输出|Text output')
@click.option('--log-file',
              help='日志文件|Log file')
def count(fastq_dir, samples, output_dir, kmer_size, threads,
          batch_size, count_separate_strands, text_output, log_file):
    """
    k-mer计数|K-mer Counting

    示例|Example: biopytools kmeria count -i fastq/ --samples samples.txt -o counts/
    """
    from biopytools.kmeria.main import main as kmeria_main

    args = ['kmeria.py', 'count',
            '-i', fastq_dir,
            '--samples', samples,
            '-o', output_dir,
            '-k', str(kmer_size),
            '-t', str(threads),
            '-b', str(batch_size)]

    if count_separate_strands:
        args.append('-C')
    if text_output:
        args.append('-T')
    if log_file:
        args.extend(['--log-file', log_file])

    original_argv = sys.argv
    sys.argv = args

    try:
        kmeria_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv


@kmeria.command('kctm')
@click.option('-i', '--input-dir', required=True,
              help='输入目录|Input directory')
@click.option('-o', '--output-dir', default='./02_kmer_matrices',
              help='输出目录|Output directory')
@click.option('-t', '--threads', default=12,
              help='线程数|Thread count (default: 12)')
@click.option('--log-file',
              help='日志文件|Log file')
def kctm(input_dir, output_dir, threads, log_file):
    """
    k-mer矩阵构建|K-mer Matrix Construction

    示例|Example: biopytools kmeria kctm -i counts/ -o matrices/
    """
    from biopytools.kmeria.main import main as kmeria_main

    args = ['kmeria.py', 'kctm',
            '-i', input_dir,
            '-o', output_dir,
            '-t', str(threads)]

    if log_file:
        args.extend(['--log-file', log_file])

    original_argv = sys.argv
    sys.argv = args

    try:
        kmeria_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv


@kmeria.command('filter')
@click.option('-i', '--input-dir', required=True,
              help='输入目录|Input directory')
@click.option('-o', '--output-dir', default='./03_filtered_matrices',
              help='输出目录|Output directory')
@click.option('-d', '--depth-file', required=True,
              help='测序深度文件|Sequencing depth file')
@click.option('-c', '--max-abund', default=1000,
              help='最大丰度|Maximum abundance (default: 1000)')
@click.option('-s', '--missing-ratio', default=0.05,
              help='缺失率|Missing ratio (default: 0.05)')
@click.option('-p', '--ploidy', default=2,
              help='倍性|Ploidy (default: 2)')
@click.option('-t', '--threads', default=12,
              help='线程数|Thread count (default: 12)')
@click.option('--log-file',
              help='日志文件|Log file')
def filter_cmd(input_dir, output_dir, depth_file, max_abund,
              missing_ratio, ploidy, threads, log_file):
    """
    k-mer矩阵过滤|K-mer Matrix Filtering

    示例|Example: biopytools kmeria filter -i matrices/ -d depth.txt -o filtered/
    """
    from biopytools.kmeria.main import main as kmeria_main

    args = ['kmeria.py', 'filter',
            '-i', input_dir,
            '-o', output_dir,
            '-d', depth_file,
            '-c', str(max_abund),
            '-s', str(missing_ratio),
            '-p', str(ploidy),
            '-t', str(threads)]

    if log_file:
        args.extend(['--log-file', log_file])

    original_argv = sys.argv
    sys.argv = args

    try:
        kmeria_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv


@kmeria.command('m2b')
@click.option('--in', '-i', 'input_dir', required=True,
              help='输入目录|Input directory')
@click.option('--out', '-o', 'output_dir', default='./04_bimbam',
              help='输出目录|Output directory')
@click.option('-t', '--threads', default=12,
              help='线程数|Thread count (default: 12)')
@click.option('--no-normalize', is_flag=True,
              help='不归一化|No normalization')
@click.option('--quantile-norm', is_flag=True,
              help='分位数归一化|Quantile normalization')
@click.option('--log-file',
              help='日志文件|Log file')
def m2b(input_dir, output_dir, threads, no_normalize, quantile_norm, log_file):
    """
    转换为BIMBAM格式|Convert to BIMBAM Format

    示例|Example: biopytools kmeria m2b -i filtered/ -o bimbam/
    """
    from biopytools.kmeria.main import main as kmeria_main

    args = ['kmeria.py', 'm2b',
            '--in', input_dir,
            '--out', output_dir,
            '-t', str(threads)]

    if no_normalize:
        args.append('--no-normalize')
    if quantile_norm:
        args.append('--quantile-norm')
    if log_file:
        args.extend(['--log-file', log_file])

    original_argv = sys.argv
    sys.argv = args

    try:
        kmeria_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv


@kmeria.command('asso')
@click.option('-i', '--input-dir', required=True,
              help='输入目录|Input directory')
@click.option('-p', '--pheno-file', required=True,
              help='表型文件|Phenotype file')
@click.option('-o', '--output-dir', default='./05_association',
              help='输出目录|Output directory')
@click.option('-n', '--pheno-col', default=1,
              help='表型列|Phenotype column (default: 1)')
@click.option('-c', '--covar-file',
              help='协变量文件|Covariate file')
@click.option('-k', '--kinship-file',
              help='亲缘关系矩阵|Kinship matrix file')
@click.option('-t', '--threads', default=12,
              help='线程数|Thread count (default: 12)')
@click.option('--log-file',
              help='日志文件|Log file')
def asso(input_dir, pheno_file, output_dir, pheno_col,
        covar_file, kinship_file, threads, log_file):
    """
    k-mer关联分析|K-mer Association Analysis

    示例|Example: biopytools kmeria asso -i bimbam/ -p pheno.txt -o results/
    """
    from biopytools.kmeria.main import main as kmeria_main

    args = ['kmeria.py', 'asso',
            '-i', input_dir,
            '-p', pheno_file,
            '-o', output_dir,
            '-n', str(pheno_col),
            '-t', str(threads)]

    if covar_file:
        args.extend(['-c', covar_file])
    if kinship_file:
        args.extend(['-k', kinship_file])
    if log_file:
        args.extend(['--log-file', log_file])

    original_argv = sys.argv
    sys.argv = args

    try:
        kmeria_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
