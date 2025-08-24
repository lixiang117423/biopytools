"""
转录组预测分析命令 | Transcriptome Prediction Analysis Command
"""

import click
import sys
from ...transcriptome_prediction.main import main as transcriptome_main


@click.command(short_help = "转录组预测分析工具")
@click.option('--genome', '-g',
              required=True,
              type=click.Path(exists=True),
              help='基因组FASTA文件 | Genome FASTA file')
@click.option('--rna-seq', '-r',
              multiple=True,
              type=click.Path(exists=True),
              help='RNA-seq FASTQ文件 (支持单端或配对末端) | RNA-seq FASTQ files (supports single-end or paired-end)')
@click.option('--samples-file',
              type=click.Path(exists=True),
              help='Trinity样本文件格式 | Trinity samples file format')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录 | Output directory')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='线程数 (默认: 88) | Number of threads (default: 88)')
# HISAT2参数
@click.option('--hisat2-min-intron',
              type=int,
              default=20,
              help='HISAT2最小内含子长度 (默认: 20) | HISAT2 minimum intron length (default: 20)')
@click.option('--hisat2-max-intron',
              type=int,
              default=500000,
              help='HISAT2最大内含子长度 (默认: 500000) | HISAT2 maximum intron length (default: 500000)')
@click.option('--hisat2-novel-splicesite',
              is_flag=True,
              help='HISAT2输出新剪接位点 | HISAT2 output novel splice sites')
@click.option('--no-dta',
              is_flag=True,
              help='禁用HISAT2的--dta选项 | Disable HISAT2 --dta option')
# StringTie参数
@click.option('--stringtie-min-length',
              type=int,
              default=200,
              help='StringTie最小转录本长度 (默认: 200) | StringTie minimum transcript length (default: 200)')
@click.option('--stringtie-min-coverage',
              type=float,
              default=1.0,
              help='StringTie最小覆盖度 (默认: 1.0) | StringTie minimum coverage (default: 1.0)')
@click.option('--stringtie-min-fpkm',
              type=float,
              default=1.0,
              help='StringTie最小FPKM (默认: 1.0) | StringTie minimum FPKM (default: 1.0)')
@click.option('--stringtie-min-iso', 
              type=float, 
              default=0.01,
              help='StringTie最小isoform (默认: 0.01) | StringTie isoform fraction (default: 0.01)')
@click.option('--stringtie-conservative',
              is_flag=True,
              help='StringTie保守组装模式 | StringTie conservative assembly mode')
# Trinity参数
@click.option('--trinity-min-contig-length',
              type=int,
              default=200,
              help='Trinity最小contig长度 (默认: 200) | Trinity minimum contig length (default: 200)')
@click.option('--trinity-max-memory',
              default='20G',
              help='Trinity最大内存使用 (默认: 20G) | Trinity maximum memory usage (default: 20G)')
@click.option('--trinity-cpu',
              type=int,
              default=88,
              help='Trinity CPU数量 (默认: 88) | Trinity CPU count (default: 88)')
@click.option('--trinity-ss-lib-type',
              type=click.Choice(['FR', 'RF', 'F', 'R']),
              help='Trinity链特异性类型 | Trinity strand-specific library type')
# PASA参数
@click.option('--pasa-max-intron-length',
              type=int,
              default=100000,
              help='PASA最大内含子长度 (默认: 100000) | PASA maximum intron length (default: 100000)')
@click.option('--pasa-min-percent-aligned',
              type=int,
              default=90,
              help='PASA最小比对百分比 (默认: 90) | PASA minimum percent aligned (default: 90)')
@click.option('--pasa-min-avg-per-id',
              type=int,
              default=95,
              help='PASA最小平均身份百分比 (默认: 95) | PASA minimum average percent identity (default: 95)')
@click.option('--pasa-aligners',
              default='gmap,blat',
              help='PASA比对器 (默认: gmap,blat) | PASA aligners (default: gmap,blat)')
@click.option('--pasa-cpu',
              type=int,
              default=88,
              help='PASA CPU数量 (默认: 88) | PASA CPU count (default: 88)')
# TransDecoder参数
@click.option('--transdecoder-min-protein-len',
              type=int,
              default=100,
              help='TransDecoder最小蛋白质长度 (默认: 100) | TransDecoder minimum protein length (default: 100)')
@click.option('--transdecoder-genetic-code',
              default='universal',
              help='TransDecoder遗传密码 (默认: universal) | TransDecoder genetic code (default: universal)')
@click.option('--transdecoder-complete-orfs-only',
              is_flag=True,
              help='TransDecoder只保留完整ORF | TransDecoder keep only complete ORFs')
# 工具路径
@click.option('--hisat2-path',
              default='hisat2',
              help='HISAT2可执行文件路径 (默认: hisat2) | HISAT2 executable path (default: hisat2)')
@click.option('--stringtie-path',
              default='stringtie',
              help='StringTie可执行文件路径 (默认: stringtie) | StringTie executable path (default: stringtie)')
@click.option('--trinity-path',
              default='Trinity',
              help='Trinity可执行文件路径 (默认: Trinity) | Trinity executable path (default: Trinity)')
@click.option('--pasa-path',
              default='Launch_PASA_pipeline.pl',
              help='PASA可执行文件路径 (默认: Launch_PASA_pipeline.pl) | PASA executable path (default: Launch_PASA_pipeline.pl)')
@click.option('--transdecoder-longorfs-path',
              default='TransDecoder.LongOrfs',
              help='TransDecoder.LongOrfs可执行文件路径 | TransDecoder.LongOrfs executable path')
@click.option('--transdecoder-predict-path',
              default='TransDecoder.Predict',
              help='TransDecoder.Predict可执行文件路径 | TransDecoder.Predict executable path')
@click.option('--samtools-path',
              default='samtools',
              help='SAMtools可执行文件路径 (默认: samtools) | SAMtools executable path (default: samtools)')
def transcriptome_prediction(**kwargs):
    """
    转录组预测分析工具
    
    集成HISAT2、StringTie、Trinity、PASA和TransDecoder的完整转录组分析流程，
    从RNA-seq数据进行转录本组装、注释和编码区预测。
    
    示例 | Examples:
    
    \b
    # 基本分析
    biopytools mrna-prediction -g genome.fa -r sample_R1.fq sample_R2.fq -o results
    
    \b
    # 使用样本文件
    biopytools mrna-prediction -g genome.fa --samples-file samples.txt -o results
    
    \b
    # 高性能自定义参数
    biopytools mrna-prediction -g genome.fa -r R1.fq R2.fq -o results \\
        -t 96 --trinity-max-memory 100G --stringtie-min-length 300
    """
    
    # 验证互斥参数
    if not kwargs.get('rna_seq') and not kwargs.get('samples_file'):
        raise click.UsageError("必须提供 --rna-seq 或 --samples-file 中的一个")
    if kwargs.get('rna_seq') and kwargs.get('samples_file'):
        raise click.UsageError("--rna-seq 和 --samples-file 不能同时使用")
    
    # 构建参数列表
    # args = ['biopytools', 'transcriptome-prediction']
    args = ['transcriptome_prediction.py']  # 改为模拟脚本名
    
    # 必需参数
    args.extend(['-g', kwargs['genome']])
    args.extend(['-o', kwargs['output']])
    
    # 输入文件
    if kwargs.get('rna_seq'):
        args.extend(['-r'] + list(kwargs['rna_seq']))
    if kwargs.get('samples_file'):
        args.extend(['--samples-file', kwargs['samples_file']])
    
    # 只在非默认值时添加参数
    defaults = {
        'threads': 88,
        'hisat2_min_intron': 20,
        'hisat2_max_intron': 500000,
        'stringtie_min_length': 200,
        'stringtie_min_coverage': 1.0,
        'stringtie_min_fpkm': 1.0,
        'trinity_min_contig_length': 200,
        'trinity_max_memory': '20G',
        'trinity_cpu': 88,
        'pasa_max_intron_length': 100000,
        'pasa_min_percent_aligned': 90,
        'pasa_min_avg_per_id': 95,
        'pasa_aligners': 'gmap,blat',
        'pasa_cpu': 88,
        'transdecoder_min_protein_len': 100,
        'transdecoder_genetic_code': 'universal',
        'hisat2_path': 'hisat2',
        'stringtie_path': 'stringtie',
        'trinity_path': 'Trinity',
        'pasa_path': 'Launch_PASA_pipeline.pl',
        'transdecoder_longorfs_path': 'TransDecoder.LongOrfs',
        'transdecoder_predict_path': 'TransDecoder.Predict',
        'samtools_path': 'samtools'
    }
    
    # 添加非默认值参数
    for key, value in kwargs.items():
        if key in ['genome', 'output', 'rna_seq', 'samples_file'] or value is None:
            continue
            
        param_name = '--' + key.replace('_', '-')
        
        if isinstance(value, bool):
            if value:
                args.append(param_name)
        elif key in defaults and value != defaults[key]:
            args.extend([param_name, str(value)])
        elif key not in defaults:
            args.extend([param_name, str(value)])
    
    # 调用原始main函数
    original_argv = sys.argv
    sys.argv = args
    
    try:
        transcriptome_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv