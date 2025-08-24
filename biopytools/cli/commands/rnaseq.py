"""
RNA-seq分析命令 | RNA-seq Analysis Command
"""

import click
import sys
from ...rnaseq.main import main as rnaseq_main


@click.command(short_help='RNA-seq表达定量分析流程',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--genome', '-g',
              required=True,
              type=click.Path(exists=True),
              help='基因组fasta文件路径 | Genome fasta file path')
@click.option('--gtf', '-f',
              required=True,
              type=click.Path(exists=True),
              help='基因注释GTF文件路径 | Gene annotation GTF file path')
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='输入fastq文件目录或样本信息文件 | Input fastq file directory or sample information file')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录 | Output directory')
@click.option('--pattern', '-p',
              type=str,
              help='Fastq文件命名模式，例如 "*.R1.fastq.gz" 或 "*_1.fq.gz"，*代表样本名 | '
                   'Fastq file naming pattern, e.g., "*.R1.fastq.gz" or "*_1.fq.gz", * represents sample name')
@click.option('--remove', '-r',
              default='no',
              type=click.Choice(['yes', 'y', 'no', 'n']),
              help='处理后删除BAM文件 | Remove BAM files after processing (default: no)')
@click.option('--threads', '-t',
              default=8,
              type=int,
              help='线程数 | Number of threads (default: 8)')
def rnaseq(genome, gtf, input, output, pattern, remove, threads):
    """
    RNA-seq分析流程：HISAT2 + StringTie (模块化版本)
    
    使用HISAT2进行转录组比对，StringTie进行转录本组装和定量，
    生成基因表达矩阵（FPKM和TPM值），适用于差异表达分析。
    
    功能特点 | Features:
    - HISAT2高效转录组比对
    - StringTie转录本组装和定量
    - 自动样本识别和批处理
    - FPKM和TPM表达量计算
    - 表达矩阵合并和标准化
    - 可选的存储空间管理
    
    输出文件 | Output Files:
    - all.fpkm.tpm.txt: 所有样本的FPKM和TPM表达矩阵
    - stringtie_output/: StringTie定量结果目录
    - fpkm_output/: 每个样本的FPKM值文件
    - *.sorted.bam: 比对结果文件（可选删除）
    - hisat2_index/: HISAT2索引文件
    
    示例 | Examples:
    
    \b
    # 基本分析
    biopytools rnaseq -g genome.fa -f genes.gtf \\
        -i /data/fastq/ -o rnaseq_results
    
    \b
    # 指定文件命名模式
    biopytools rnaseq --genome reference.fasta --gtf annotation.gtf \\
        --input ./samples/ --output ./analysis/ \\
        --pattern "*.R1.fastq.gz"
    
    \b
    # 高线程数分析并删除BAM文件
    biopytools rnaseq -g genome.fa -f genes.gtf \\
        -i /data/rna_samples/ -o results/ \\
        -t 32 --remove yes
    
    \b
    # 自定义模式的配对端测序数据
    biopytools rnaseq -g reference.fasta -f transcripts.gtf \\
        -i /project/fastq_files/ -o /project/results/ \\
        --pattern "*_1.clean.fq.gz" --threads 16
    
    \b
    # 保留BAM文件用于后续分析
    biopytools rnaseq -g genome.fasta -f annotation.gtf \\
        -i ./raw_data/ -o ./expression_analysis/ \\
        --remove no --threads 24
    
    文件命名模式说明 | File Naming Pattern:
    
    使用通配符*表示样本名称，工具会自动识别配对的R1和R2文件：
    - "*.R1.fastq.gz" 匹配: sample1.R1.fastq.gz, sample1.R2.fastq.gz
    - "*_1.fq.gz" 匹配: sample_1.fq.gz, sample_2.fq.gz  
    - "*_R1.clean.fastq" 匹配: sample_R1.clean.fastq, sample_R2.clean.fastq
    
    如果不指定pattern，工具会尝试自动检测常见的命名模式。
    
    输入目录结构示例 | Input Directory Structure:
    /data/fastq/
    ├── sample1_R1.fastq.gz
    ├── sample1_R2.fastq.gz
    ├── sample2_R1.fastq.gz
    ├── sample2_R2.fastq.gz
    └── ...
    
    输出目录结构 | Output Directory Structure:
    results/
    ├── all.fpkm.tpm.txt           # 主要结果文件
    ├── hisat2_index/              # HISAT2索引
    ├── stringtie_output/          # StringTie结果
    │   ├── sample1.gtf
    │   └── sample2.gtf
    ├── fpkm_output/               # FPKM值文件
    │   ├── sample1.fpkm.txt
    │   └── sample2.fpkm.txt
    └── *.sorted.bam               # BAM文件（可选）
    
    分析流程 | Analysis Pipeline:
    1. 构建HISAT2基因组索引
    2. 解析输入样本和文件配对
    3. 对每个样本执行：
       a. HISAT2转录组比对
       b. StringTie转录本定量
       c. 提取基因表达值（FPKM）
       d. 可选删除BAM文件节省空间
    4. 合并所有样本的表达矩阵
    5. 计算TPM值并标准化
    6. 生成分析总结报告
    
    表达量指标说明 | Expression Metrics:
    - FPKM: Fragments Per Kilobase per Million，考虑基因长度和测序深度
    - TPM: Transcripts Per Million，标准化后便于样本间比较
    
    适用场景 | Use Cases:
    - 差异基因表达分析
    - 转录组比较研究
    - 基因表达谱构建
    - 功能富集分析前处理
    - 共表达网络分析
    
    性能建议 | Performance Tips:
    - 增加线程数可显著提升比对速度
    - SSD存储有助于提高I/O性能
    - 大项目可选择删除BAM文件节省空间
    - 确保有足够内存用于索引构建
    - 基因组索引可重复使用，避免重复构建
    
    质量控制建议 | Quality Control:
    - 确保FASTQ文件质量良好
    - 检查基因组和GTF文件版本一致性
    - 验证样本命名的一致性
    - 监控比对率和定量质量
    """
    
    # 构建参数列表传递给原始main函数
    args = ['rnaseq.py']
    
    # 必需参数
    args.extend(['-g', genome])
    args.extend(['-f', gtf])
    args.extend(['-i', input])
    args.extend(['-o', output])
    
    # 可选参数（只在非默认值时添加）
    if pattern:
        args.extend(['-p', pattern])
    
    if remove != 'no':
        args.extend(['-r', remove])
    
    if threads != 8:
        args.extend(['-t', str(threads)])
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        rnaseq_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n分析被用户中断 | Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"分析失败 | Analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv