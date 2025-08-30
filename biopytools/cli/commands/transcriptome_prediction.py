"""
🧬 转录组预测分析命令 | Transcriptome Prediction Analysis Command
"""

import click
import sys
from ...transcriptome_prediction.main import main as transcriptome_main


@click.command(short_help='🧬 转录组预测分析工具：集成RNA-seq转录本组装、注释和编码区预测',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--genome', '-g',
              required=True,
              type=click.Path(exists=True),
              help='🧬 基因组FASTA文件 | Genome FASTA file')
@click.option('--rna-seq', '-r',
              multiple=True,
              type=click.Path(exists=True),
              help='📊 RNA-seq FASTQ文件 (支持单端或配对末端) | RNA-seq FASTQ files (supports single-end or paired-end)')
@click.option('--samples-file',
              type=click.Path(exists=True),
              help='📋 Trinity样本文件格式 | Trinity samples file format')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='📁 输出目录 | Output directory')
# 流程控制参数
# @click.option('--resume',
#               is_flag=True,
#               default=True,
#               help='🔄 启用断点续传，跳过已完成的步骤（默认启用） | Enable resume mode, skip completed steps (enabled by default)')
@click.option('--no-resume',
              is_flag=True,
              help='🚫 禁用断点续传，重新运行所有步骤 | Disable resume mode, rerun all steps')
@click.option('--skip-trinity',
              is_flag=True,
              help='🔗 跳过Trinity de novo组装步骤 | Skip Trinity de novo assembly step')
@click.option('--step',
              type=click.Choice(['alignment', 'stringtie', 'trinity', 'pasa', 'transdecoder']),
              help='🎯 只运行特定步骤 | Run only specific step')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🧵 线程数 (默认: 88) | Number of threads (default: 88)')
# HISAT2参数
@click.option('--hisat2-min-intron',
              type=int,
              default=20,
              help='🎯 HISAT2最小内含子长度 (默认: 20) | HISAT2 minimum intron length (default: 20)')
@click.option('--hisat2-max-intron',
              type=int,
              default=500000,
              help='🎯 HISAT2最大内含子长度 (默认: 500000) | HISAT2 maximum intron length (default: 500000)')
@click.option('--hisat2-novel-splicesite',
              is_flag=True,
              help='🎯 HISAT2输出新剪接位点 | HISAT2 output novel splice sites')
@click.option('--no-dta',
              is_flag=True,
              help='🎯 禁用HISAT2的--dta选项 | Disable HISAT2 --dta option')
# StringTie参数
@click.option('--stringtie-min-length',
              type=int,
              default=200,
              help='🧩 StringTie最小转录本长度 (默认: 200) | StringTie minimum transcript length (default: 200)')
@click.option('--stringtie-min-coverage',
              type=float,
              default=1.0,
              help='🧩 StringTie最小覆盖度 (默认: 1.0) | StringTie minimum coverage (default: 1.0)')
@click.option('--stringtie-min-fpkm',
              type=float,
              default=1.0,
              help='🧩 StringTie最小FPKM (默认: 1.0) | StringTie minimum FPKM (default: 1.0)')
@click.option('--stringtie-min-iso',
              type=float,
              default=0.01,
              help='🧩 StringTie最小isoform比例 (默认: 0.01) | StringTie isoform fraction (default: 0.01)')
@click.option('--stringtie-conservative',
              is_flag=True,
              help='🧩 StringTie保守组装模式 | StringTie conservative assembly mode')
# Trinity参数
@click.option('--trinity-min-contig-length',
              type=int,
              default=200,
              help='🔗 Trinity最小contig长度 (默认: 200) | Trinity minimum contig length (default: 200)')
@click.option('--trinity-max-memory',
              default='200G',
              help='🔗 Trinity最大内存使用 (默认: 200G) | Trinity maximum memory usage (default: 200G)')
@click.option('--trinity-cpu',
              type=int,
              default=88,
              help='🔗 Trinity CPU数量 (默认: 88) | Trinity CPU count (default: 88)')
@click.option('--trinity-ss-lib-type',
              type=click.Choice(['FR', 'RF', 'F', 'R']),
              help='🔗 Trinity链特异性类型 | Trinity strand-specific library type')
# PASA参数
@click.option('--pasa-max-intron-length',
              type=int,
              default=100000,
              help='📍 PASA最大内含子长度 (默认: 100000) | PASA maximum intron length (default: 100000)')
@click.option('--pasa-min-percent-aligned',
              type=int,
              default=90,
              help='📍 PASA最小比对百分比 (默认: 90) | PASA minimum percent aligned (default: 90)')
@click.option('--pasa-min-avg-per-id',
              type=int,
              default=95,
              help='📍 PASA最小平均身份百分比 (默认: 95) | PASA minimum average percent identity (default: 95)')
@click.option('--pasa-aligners',
              default='gmap,blat',
              help='📍 PASA比对器 (默认: gmap,blat) | PASA aligners (default: gmap,blat)')
@click.option('--pasa-cpu',
              type=int,
              default=88,
              help='📍 PASA CPU数量 (默认: 88) | PASA CPU count (default: 88)')
# TransDecoder参数
@click.option('--transdecoder-min-protein-len',
              type=int,
              default=100,
              help='🔍 TransDecoder最小蛋白质长度 (默认: 100) | TransDecoder minimum protein length (default: 100)')
@click.option('--transdecoder-genetic-code',
              default='universal',
              help='🔍 TransDecoder遗传密码 (默认: universal) | TransDecoder genetic code (default: universal)')
@click.option('--transdecoder-complete-orfs-only',
              is_flag=True,
              help='🔍 TransDecoder只保留完整ORF | TransDecoder keep only complete ORFs')
# 工具路径
@click.option('--hisat2-path',
              default='hisat2',
              help='🛠️ HISAT2可执行文件路径 (默认: hisat2) | HISAT2 executable path (default: hisat2)')
@click.option('--stringtie-path',
              default='stringtie',
              help='🛠️ StringTie可执行文件路径 (默认: stringtie) | StringTie executable path (default: stringtie)')
@click.option('--trinity-path',
              default='Trinity',
              help='🛠️ Trinity可执行文件路径 (默认: Trinity) | Trinity executable path (default: Trinity)')
@click.option('--pasa-path',
              default='Launch_PASA_pipeline.pl',
              help='🛠️ PASA可执行文件路径 (默认: Launch_PASA_pipeline.pl) | PASA executable path (default: Launch_PASA_pipeline.pl)')
@click.option('--transdecoder-longorfs-path',
              default='TransDecoder.LongOrfs',
              help='🛠️ TransDecoder.LongOrfs可执行文件路径 | TransDecoder.LongOrfs executable path')
@click.option('--transdecoder-predict-path',
              default='TransDecoder.Predict',
              help='🛠️ TransDecoder.Predict可执行文件路径 | TransDecoder.Predict executable path')
@click.option('--samtools-path',
              default='samtools',
              help='🛠️ SAMtools可执行文件路径 (默认: samtools) | SAMtools executable path (default: samtools)')
def transcriptome_prediction(**kwargs):
    """
    🧬 转录组预测分析工具 | Transcriptome-based prediction analysis tool
    
    一个完整的转录组分析流程工具，集成HISAT2、StringTie、Trinity、PASA和TransDecoder，
    从RNA-seq数据进行转录本组装、注释和编码区预测，适用于转录组学研究中的基因发现和功能注释。
    
    功能特点 | Features:
    - 🎯 全自动HISAT2基因组比对和转录本组装
    - 🧩 StringTie参考基因组导向转录本重构
    - 🔗 Trinity de novo转录本组装（可选）
    - 📍 PASA基因结构综合注释和验证
    - 🔍 TransDecoder蛋白质编码区预测
    - 🔄 智能断点续传，支持大规模数据处理
    - 📊 详细的分析日志记录和质量统计
    - ⚙️ 模块化可扩展架构设计
    
    分析流程 | Analysis Pipeline:
    1. 🎯 HISAT2比对：将RNA-seq读段比对到参考基因组
    2. 🧩 StringTie组装：重构转录本并合并多样本结果
    3. 🔗 Trinity组装：de novo转录本组装（可选）
    4. 📍 PASA注释：整合转录证据进行基因结构注释
    5. 🔍 TransDecoder预测：识别和预测蛋白质编码区
    
    应用场景 | Use Cases:
    - 新物种转录组注释和基因发现
    - 转录本变异体识别和功能分析
    - 非模式生物基因组功能注释
    - 转录组辅助基因组注释改进
    - 比较转录组学和进化分析
    - 个性化医学中的转录组分析
    
    示例 | Examples:
    
    \b
    # 基本转录组分析
    biopytools transcriptome-prediction -g genome.fa -r sample_R1.fq sample_R2.fq -o results
    
    \b
    # 跳过Trinity步骤的快速分析
    biopytools transcriptome-prediction -g genome.fa -r R1.fq R2.fq -o results --skip-trinity
    
    \b
    # 使用样本文件进行多样本分析
    biopytools transcriptome-prediction -g genome.fa --samples-file samples.txt -o results
    
    \b
    # 禁用断点续传，重新运行所有步骤
    biopytools transcriptome-prediction -g genome.fa -r R1.fq R2.fq -o results --no-resume
    
    \b
    # 只运行特定步骤
    biopytools transcriptome-prediction -g genome.fa -r R1.fq R2.fq -o results --step alignment
    biopytools transcriptome-prediction -g genome.fa -r R1.fq R2.fq -o results --step stringtie
    biopytools transcriptome-prediction -g genome.fa -r R1.fq R2.fq -o results --step trinity
    
    \b
    # 高性能服务器优化配置
    biopytools transcriptome-prediction -g genome.fa -r R1.fq R2.fq -o results \\
        -t 128 --trinity-max-memory 500G --trinity-cpu 128 \\
        --pasa-cpu 128 --stringtie-min-length 300
    
    \b
    # 自定义质量参数的严格分析
    biopytools transcriptome-prediction -g genome.fa -r R1.fq R2.fq -o results \\
        --stringtie-min-coverage 2.0 --stringtie-min-fpkm 2.0 \\
        --transdecoder-min-protein-len 150 --stringtie-conservative
    
    \b
    # 指定外部工具路径
    biopytools transcriptome-prediction -g genome.fa -r R1.fq R2.fq -o results \\
        --hisat2-path /usr/local/bin/hisat2 \\
        --trinity-path /opt/Trinity/Trinity \\
        --pasa-path /opt/PASApipeline/Launch_PASA_pipeline.pl
    
    \b
    # 完整高级分析流程示例
    biopytools transcriptome-prediction -g reference_genome.fa \\
        -r tissue1_R1.fq tissue1_R2.fq tissue2_R1.fq tissue2_R2.fq \\
        -o comprehensive_transcriptome \\
        -t 96 --trinity-max-memory 300G \\
        --hisat2-max-intron 1000000 \\
        --stringtie-min-length 250 --stringtie-min-coverage 1.5 \\
        --pasa-min-percent-aligned 95 \\
        --transdecoder-min-protein-len 120 \\
        --resume
    
    输入文件格式 | Input File Formats:
    
    基因组文件要求:
    - 标准FASTA格式
    - 已完成组装的基因组序列
    - 建议包含完整的染色体序列
    - 支持压缩格式(.fa.gz)
    
    RNA-seq文件格式:
    - 标准FASTQ格式 (支持.fq, .fastq, .fq.gz, .fastq.gz)
    - 支持单端测序 (Single-end) 数据
    - 支持双端测序 (Paired-end) 数据
    - 配对文件需按R1, R2顺序提供
    - 支持多样本同时分析
    
    样本文件格式 (Trinity格式):
    - 制表符分隔的文本文件
    - 每行格式: condition_name[TAB]replicate_name[TAB]path_to_left_reads[TAB]path_to_right_reads
    - 单端数据省略右端文件列
    - 支持生物学重复和技术重复
    
    样本文件示例:
    control    rep1    ctrl_1_R1.fq    ctrl_1_R2.fq
    control    rep2    ctrl_2_R1.fq    ctrl_2_R2.fq
    treatment  rep1    treat_1_R1.fq   treat_1_R2.fq
    treatment  rep2    treat_2_R1.fq   treat_2_R2.fq
    
    工具参数说明 | Tool Parameters:
    
    HISAT2比对参数:
    --hisat2-min-intron: 最小内含子长度，影响剪接位点识别
    --hisat2-max-intron: 最大内含子长度，根据物种调整
    --hisat2-novel-splicesite: 输出新发现的剪接位点
    --no-dta: 禁用downstream transcriptome analysis模式
    
    StringTie组装参数:
    --stringtie-min-length: 转录本最小长度阈值
    --stringtie-min-coverage: 最小reads覆盖度
    --stringtie-min-fpkm: 最小FPKM表达阈值
    --stringtie-min-iso: isoform最小比例阈值
    --stringtie-conservative: 保守模式，减少假阳性
    
    Trinity参数:
    --trinity-min-contig-length: 输出contig最小长度
    --trinity-max-memory: 最大内存使用限制
    --trinity-cpu: 使用的CPU核心数
    --trinity-ss-lib-type: 链特异性文库类型
    
    PASA注释参数:
    --pasa-max-intron-length: 允许的最大内含子长度
    --pasa-min-percent-aligned: 最小比对覆盖百分比
    --pasa-min-avg-per-id: 最小平均序列相似性
    --pasa-aligners: 使用的比对工具组合
    
    TransDecoder参数:
    --transdecoder-min-protein-len: 预测蛋白质最小长度
    --transdecoder-genetic-code: 使用的遗传密码表
    --transdecoder-complete-orfs-only: 只保留完整的开放阅读框
    
    输出文件说明 | Output Files:
    
    核心结果文件:
    - hisat2_alignment/: HISAT2比对结果BAM文件
    - stringtie_assembly/: StringTie转录本GTF文件
    - trinity_assembly/: Trinity de novo组装FASTA文件
    - pasa_annotation/: PASA注释GFF3/GTF文件
    - transdecoder_prediction/: TransDecoder蛋白质序列和注释
    
    详细输出结构:
    output_directory/
    ├── hisat2_alignment/
    │   ├── sample1_sorted.bam
    │   ├── sample2_sorted.bam
    │   └── alignment_stats.txt
    ├── stringtie_assembly/
    │   ├── sample1.gtf
    │   ├── sample2.gtf
    │   └── merged.gtf
    ├── trinity_assembly/
    │   └── Trinity.fasta
    ├── pasa_annotation/
    │   ├── pasa_assemblies.gff3
    │   └── pasa_assemblies.gtf
    ├── transdecoder_prediction/
    │   ├── transcripts.fasta.transdecoder.cds
    │   ├── transcripts.fasta.transdecoder.pep
    │   └── transcripts.fasta.transdecoder.gff3
    └── transcriptome_analysis_report.txt
    
    断点续传功能 | Resume Functionality:
    
    系统自动检测已完成的步骤:
    - ✅ 检查BAM文件存在性（比对步骤）
    - ✅ 检查GTF文件完整性（StringTie步骤）
    - ✅ 检查Trinity组装文件（Trinity步骤）
    - ✅ 检查PASA注释结果（PASA步骤）
    - ✅ 检查TransDecoder预测文件（预测步骤）
    
    断点续传优势:
    - 💾 节省计算资源和时间
    - 🔄 支持大规模数据分析中断恢复
    - 📊 保持数据一致性和完整性
    - ⚡ 快速跳过耗时的组装步骤
    
    性能和系统要求 | Performance & System Requirements:
    
    依赖软件:
    - HISAT2 (v2.1.0+): 快速RNA-seq比对工具
    - StringTie (v2.0+): 转录本组装和定量
    - Trinity (v2.8.0+): de novo转录组组装
    - PASA (v2.4.0+): 转录本注释流程
    - TransDecoder (v5.0+): 编码区预测工具
    - SAMtools (v1.9+): 序列数据处理
    
    系统建议:
    - RAM: 至少32GB，推荐128GB+用于大基因组
    - CPU: 多核处理器，推荐64核+用于Trinity
    - 存储: 至少输入数据大小的10倍自由空间
    - 网络: 高速存储I/O用于大文件处理
    
    数据规模估算:
    - 小型基因组(100Mb)，10M reads: ~8GB内存，~50GB存储
    - 中型基因组(1Gb)，50M reads: ~64GB内存，~500GB存储
    - 大型基因组(3Gb)，200M reads: ~256GB内存，~2TB存储
    - 复杂基因组，超大数据集: 考虑分批处理和集群计算
    
    故障排除 | Troubleshooting:
    
    常见问题:
    1. "Tool not found": 检查工具安装和PATH设置
    2. "Memory error": 增加系统内存或调整参数
    3. "Disk space": 确保足够的存储空间
    4. "Trinity failed": 检查内存限制和CPU核心数
    5. "PASA alignment failed": 验证输入文件格式和参数
    
    优化建议:
    - 根据数据量调整内存和CPU参数
    - 使用SSD存储提高I/O性能
    - 大数据集考虑预先质控和过滤
    - 监控系统资源使用情况
    - 设置合理的质量阈值参数
    
    最佳实践 | Best Practices:
    
    1. 数据准备:
       - 使用高质量的RNA-seq数据
       - 预先进行序列质控和清理
       - 确保参考基因组质量良好
       - 合理设计生物学重复实验
    
    2. 参数调优:
       - 根据物种特性调整内含子长度
       - 基于数据质量设置阈值参数
       - 考虑计算资源限制优化配置
       - 测试不同参数对结果的影响
    
    3. 结果验证:
       - 检查转录本组装质量统计
       - 验证预测基因的功能注释
       - 比较不同方法的结果一致性
       - 使用已知基因验证预测准确性
    
    4. 下游分析:
       - 功能注释和通路富集分析
       - 差异表达基因识别
       - 转录本变异体分析
       - 基因组注释文件生成
    
    引用和参考 | Citation & References:
    
    如果在学术研究中使用此工具，请引用相关的方法学文献:
    - HISAT2: Kim et al. (2019) Nature Biotechnology
    - StringTie: Pertea et al. (2015) Nature Biotechnology  
    - Trinity: Grabherr et al. (2011) Nature Biotechnology
    - PASA: Haas et al. (2003) Nucleic Acids Research
    - TransDecoder: Haas et al. (2013) Nature Protocols
    """
    
    # 验证输入参数
    if not kwargs.get('rna_seq') and not kwargs.get('samples_file'):
        raise click.UsageError("❌ 必须提供 --rna-seq 或 --samples-file 中的一个 | Must provide either --rna-seq or --samples-file")
    
    if kwargs.get('rna_seq') and kwargs.get('samples_file'):
        raise click.UsageError("❌ --rna-seq 和 --samples-file 不能同时使用 | Cannot use --rna-seq and --samples-file together")
    
    # 验证步骤和跳过Trinity的冲突
    if kwargs.get('step') == 'trinity' and kwargs.get('skip_trinity'):
        raise click.UsageError("❌ 无法运行Trinity步骤：指定了--skip-trinity参数 | Cannot run Trinity step: --skip-trinity specified")
    
    # 构建参数列表传递给原始main函数
    args = ['transcriptome_prediction.py']
    
    # 必需参数
    args.extend(['-g', kwargs['genome']])
    args.extend(['-o', kwargs['output']])
    
    # 输入文件参数
    if kwargs.get('rna_seq'):
        args.extend(['-r'] + list(kwargs['rna_seq']))
    if kwargs.get('samples_file'):
        args.extend(['--samples-file', kwargs['samples_file']])
    
    # 流程控制参数
    if kwargs.get('no_resume'):
        args.append('--no-resume')
    elif not kwargs.get('resume', True):  # 如果明确设置为False
        args.append('--no-resume')
    
    if kwargs.get('skip_trinity'):
        args.append('--skip-trinity')
    
    if kwargs.get('step'):
        args.extend(['--step', kwargs['step']])
    
    # 定义默认值，只在非默认值时添加参数
    defaults = {
        'threads': 88,
        'hisat2_min_intron': 20,
        'hisat2_max_intron': 500000,
        'stringtie_min_length': 200,
        'stringtie_min_coverage': 1.0,
        'stringtie_min_fpkm': 1.0,
        'stringtie_min_iso': 0.01,
        'trinity_min_contig_length': 200,
        'trinity_max_memory': '200G',
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
    
    # 添加其他参数
    skip_params = {'genome', 'output', 'rna_seq', 'samples_file', 'resume', 'no_resume', 'skip_trinity', 'step'}
    
    for key, value in kwargs.items():
        if key in skip_params or value is None:
            continue
            
        param_name = '--' + key.replace('_', '-')
        
        if isinstance(value, bool):
            if value:  # 只在True时添加flag
                args.append(param_name)
        elif key in defaults:
            if value != defaults[key]:  # 只在非默认值时添加
                args.extend([param_name, str(value)])
        else:
            # 不在默认值列表中的参数直接添加
            args.extend([param_name, str(value)])
    
    # 保存并替换sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始main函数
        transcriptome_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 转录组预测分析被用户中断 | Transcriptome prediction analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 转录组预测分析失败 | Transcriptome prediction analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv