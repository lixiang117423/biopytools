"""
序列提取工具命令 | Sequence Extraction Tool Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_sequence_main():
    """懒加载序列提取工具main函数 | Lazy load sequence extraction main function"""
    try:
        from ...seq_extractor.main import main as sequence_main
        return sequence_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件是否存在（仅在非帮助模式下）| Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"❌ 文件不存在 | File does not exist: {file_path}")
    return file_path


@click.command(short_help='🧬 序列提取工具：根据区域文件从FASTA文件中提取指定序列片段',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--sequence-file', '-s',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='🧬 输入序列文件 (FASTA格式) | Input sequence file (FASTA format)')
@click.option('--regions-file', '-r',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='📊 区域文件 (类似BED格式: 染色体 起始 终止) | Regions file (BED-like format: chromosome start end)')
@click.option('--output-file', '-o',
              required=True,
              type=click.Path(),
              help='💾 输出序列文件 | Output sequence file')
@click.option('--type', '--sequence-type',
              type=click.Choice(['dna', 'protein']),
              default='dna',
              help='🔬 序列类型 | Sequence type (default: dna)')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='🚀 线程数 | Number of threads (default: 88)')
@click.option('--merge-output',
              is_flag=True,
              default=True,
              help='📦 合并输出到一个文件 | Merge output to one file (default: enabled)')
@click.option('--separate-output',
              is_flag=True,
              help='📂 分别输出到多个文件 | Output to separate files')
@click.option('--no-headers',
              is_flag=True,
              help='🏷️ 不包含区域信息在序列名中 | Do not include region info in sequence names')
@click.option('--reverse-complement',
              is_flag=True,
              help='🔄 反向互补DNA序列 | Reverse complement DNA sequences')
@click.option('--translate',
              is_flag=True,
              help='🔄 将DNA翻译为蛋白质 | Translate DNA to protein')
@click.option('--line-width',
              type=int,
              default=80,
              help='📏 FASTA序列每行字符数 | Characters per line in FASTA sequence (default: 80)')
@click.option('--verbose', '-v',
              is_flag=True,
              default=True,
              help='📝 详细输出 | Verbose output (default: enabled)')
@click.option('--quiet', '-q',
              is_flag=True,
              help='🔇 安静模式 | Quiet mode')
@click.option('--samtools-path',
              default='samtools',
              help='🛠️  samtools软件路径 | samtools software path (default: samtools)')
def parse_seq(sequence_file, regions_file, output_file, type, threads, merge_output,
                      separate_output, no_headers, reverse_complement, translate, line_width,
                      verbose, quiet, samtools_path):
    """
    🧬 序列提取工具 | Sequence Extraction Tool
    
    一个专业的序列片段提取工具，用于根据指定的基因组区域从FASTA文件中精确提取
    DNA或蛋白质序列片段，支持批量处理、多线程优化和多种序列处理功能，
    广泛应用于基因组学、转录组学和蛋白质组学研究。
    
    ✨ 功能特点 | Features:
    - 🎯 精确的基因组区域序列提取
    - 🚀 高效的多线程并行处理
    - 🔄 DNA反向互补和蛋白质翻译
    - 📊 灵活的输出格式和组织方式
    - 🧬 支持DNA和蛋白质序列类型
    - 📋 批量区域处理能力
    - 🛠️ 与samtools等工具集成
    
    🔄 处理流程 | Processing Pipeline:
    1. 📖 解析输入FASTA序列文件
    2. 📊 读取基因组区域坐标文件
    3. 🎯 根据坐标提取对应序列片段
    4. 🔄 应用序列处理操作(翻译/反向互补)
    5. 🏷️ 生成带有区域信息的序列标识
    6. 💾 输出处理后的FASTA序列
    
    🎯 应用场景 | Use Cases:
    - 🧬 基因组特定区域序列获取
    - 🔍 功能基因序列提取和分析
    - 📊 基因组注释区域序列收集
    - 🧪 PCR引物设计序列准备
    - 🌳 系统发育分析序列整理
    - 🔬 比较基因组学序列对比
    - 💾 基因组数据库序列提取
    
    💡 示例 | Examples:
    
    \b
    # 🎯 基本DNA序列提取
    biopytools parse-seq -s genome.fasta -r regions.bed \\
        -o extracted.fasta --type dna
    
    \b
    # 🧪 蛋白质序列提取
    biopytools parse-seq -s proteins.fasta -r regions.bed \\
        -o extracted_proteins.fasta --type protein
    
    \b
    # 🚀 高性能多线程处理
    biopytools parse-seq -s large_genome.fasta -r regions.bed \\
        -o extracted.fasta --type dna --threads 64
    
    \b
    # 🔄 反向互补DNA序列提取
    biopytools parse-seq -s genome.fasta -r regions.bed \\
        -o reverse_comp.fasta --type dna --reverse-complement
    
    \b
    # 🔄 DNA翻译为蛋白质
    biopytools parse-seq -s cds_sequences.fasta -r regions.bed \\
        -o translated.fasta --type dna --translate
    
    \b
    # 📂 分别输出到多个文件
    biopytools parse-seq -s genome.fasta -r regions.bed \\
        -o extracted --separate-output --no-headers
    
    \b
    # 🔇 安静模式处理
    biopytools parse-seq -s genome.fasta -r regions.bed \\
        -o extracted.fasta --quiet --line-width 100
    
    📁 输入文件要求 | Input File Requirements:
    
    🧬 序列文件格式 (FASTA):
    - 📄 标准FASTA格式文件 (.fa, .fasta, .fas)
    - 🧬 DNA序列：只含ATCG字符
    - 🧪 蛋白质序列：标准20个氨基酸字符
    - 🏷️ 序列标识符应唯一且有意义
    - 📦 支持压缩文件 (.gz, .bz2)
    
    📊 区域文件格式 (BED-like):
    制表符分隔的文本文件，每行包含：
    - 🏷️ 列1：序列标识符 (对应FASTA中的序列名)
    - 📍 列2：起始位置 (0-based，包含)
    - 📍 列3：终止位置 (0-based，不包含)
    - 📝 列4：可选的区域名称或描述
    
    📄 区域文件示例:
    chr1    1000    2000    gene1
    chr1    5000    6000    gene2
    chr2    3000    4500    gene3
    mitochondria    100    500    cox1
    
    🧬 FASTA文件示例:
    >chr1
    ATCGATCGATCGATCGATCG...
    >chr2
    GCTAGCTAGCTAGCTAGCTA...
    >mitochondria
    TTAATTAATTAATTAATTAA...
    
    ⚙️ 参数详细说明 | Parameter Details:
    
    📋 必需参数:
    --sequence-file (-s): 🧬 输入FASTA序列文件
        - 包含待提取序列的完整FASTA文件
        - 序列标识符需与区域文件中的标识符对应
    
    --regions-file (-r): 📊 基因组区域坐标文件
        - BED格式或类似的坐标文件
        - 定义需要提取的序列区间
    
    --output-file (-o): 💾 输出文件路径
        - 提取序列的输出FASTA文件
        - 支持绝对路径和相对路径
    
    🔬 序列类型:
    --type: 🧬 指定序列数据类型
        - dna: DNA序列，支持反向互补和翻译
        - protein: 蛋白质序列，仅支持直接提取
    
    ⚡ 性能参数:
    --threads (-t): 🚀 并行处理线程数
        - 建议设置为CPU核心数
        - 影响大文件处理速度
        - 多区域提取时显著提升性能
    
    📂 输出控制:
    --merge-output: 📦 合并所有提取序列到单一文件
        - 默认启用，适合后续分析
        - 生成统一的FASTA输出文件
    
    --separate-output: 📂 为每个区域生成独立文件
        - 与merge-output互斥
        - 适合需要分别处理各序列的场景
    
    --no-headers: 🏷️ 简化序列标识符
        - 不在序列名中包含区域坐标信息
        - 生成更简洁的FASTA标题
    
    🔄 DNA序列处理选项:
    --reverse-complement: 🔄 生成反向互补序列
        - 仅适用于DNA序列类型
        - 获取DNA双链的另一条链序列
        - 常用于引物设计和基因分析
    
    --translate: 🔄 将DNA翻译为蛋白质
        - 使用标准遗传密码表
        - 自动处理起始和终止密码子
        - 输出对应的氨基酸序列
    
    📏 格式化选项:
    --line-width: 📏 FASTA序列行宽设置
        - 控制输出序列每行的字符数
        - 标准值为60或80字符
        - 影响文件可读性和兼容性
    
    📝 日志控制:
    --verbose (-v): 📝 详细日志输出模式
        - 默认启用，显示处理进度
        - 输出详细的操作信息
    
    --quiet (-q): 🔇 静默处理模式
        - 与verbose互斥
        - 仅输出错误和警告信息
    
    🛠️ 工具路径:
    --samtools-path: 🛠️ samtools程序路径
        - 用于处理压缩文件和索引
        - 支持BAM/CRAM格式输入(扩展功能)
    
    📁 输出文件说明 | Output Files:
    
    📦 合并输出模式:
    单一FASTA文件包含所有提取的序列：
    >region_1_chr1:1000-2000
    ATCGATCGATCGATCG...
    >region_2_chr1:5000-6000
    GCTAGCTAGCTAGCTA...
    
    📂 分离输出模式:
    为每个区域生成独立文件：
    - 📄 {output_prefix}_region1.fasta
    - 📄 {output_prefix}_region2.fasta
    - 📄 {output_prefix}_region3.fasta
    
    🏷️ 序列标识符格式:
    标准模式：region_name_seqid:start-end
    简化模式(--no-headers)：region_name
    
    ⚡ 性能和系统要求 | Performance & System Requirements:
    
    💻 系统建议:
    - 💾 RAM: 至少2GB，大基因组建议8GB+
    - 🗄️ 存储: 输出目录需要输入文件1-2倍空间
    - 🧵 CPU: 多核处理器提升并行提取效率
    - 🐍 Python: 3.7+版本
    
    📊 处理能力估算:
    - 📄 小文件(<100MB): 几秒到几分钟
    - 📋 中等文件(100MB-1GB): 几分钟到半小时
    - 📁 大文件(1GB-10GB): 半小时到几小时
    - 🚀 多线程可显著提升处理速度
    
    💾 内存使用估算:
    - 🔵 基础内存: ~200MB
    - 📈 序列相关: FASTA文件大小的1-2倍
    - 🧵 线程相关: 每线程额外~50MB
    
    🛠️ 故障排除 | Troubleshooting:
    
    ⚠️ 常见问题:
    1. ❌ "Sequence not found": 检查区域文件中的序列标识符
    2. 📍 "Invalid coordinates": 验证起始终止位置合理性
    3. 💾 "内存不足": 减少线程数或处理更小文件
    4. 📄 "文件格式错误": 确认FASTA和区域文件格式
    5. 🔧 "samtools error": 检查samtools安装和路径
    
    📋 数据质量检查:
    - 🔍 验证FASTA文件格式完整性
    - 📊 检查区域坐标范围合理性
    - 🏷️ 确认序列标识符一致性
    - 📏 验证提取序列长度正确性
    
    💡 优化建议:
    - 🎯 合理设置线程数平衡速度和内存
    - 📊 预先验证输入文件格式和完整性
    - 🗄️ 使用高速存储提升I/O性能
    - 📈 监控系统资源使用情况
    
    🏆 最佳实践 | Best Practices:
    
    1️⃣ 📊 数据准备:
       - ✅ 使用标准格式的FASTA和BED文件
       - 🔍 预先验证文件完整性和格式
       - 🏷️ 确保序列标识符命名一致性
       - 📍 检查区域坐标的生物学合理性
    
    2️⃣ ⚙️ 参数选择:
       - 🧵 根据系统配置设置合适线程数
       - 🔬 根据数据类型选择正确的序列类型
       - 📂 考虑下游分析选择输出格式
       - 🔄 DNA序列根据需要选择处理选项
    
    3️⃣ 🔍 质量控制:
       - 📏 验证提取序列长度符合预期
       - 🎲 随机检查部分序列提取准确性
       - 📊 比较提取序列数量与区域数量
       - 🧬 验证DNA翻译结果的生物学意义
    
    4️⃣ 📦 结果处理:
       - 💾 及时备份重要的提取结果
       - 📋 记录详细的处理参数和日志
       - 🗂️ 建立清晰的文件命名规范
       - 📊 保留原始数据便于验证和重现
    
    🔍 结果验证建议 | Result Validation:
    
    🤖 自动验证:
    - 📊 提取序列数量与区域数量对比
    - 📏 序列长度与预期区间长度检查
    - 🧬 DNA序列ATCG字符比例验证
    - 🏷️ 输出文件格式规范性检查
    
    👥 人工验证:
    - 🎲 随机选择区域验证提取准确性
    - 🔍 检查关键基因区域序列完整性
    - 📚 与已知序列数据库进行比对
    - 🧪 验证翻译蛋白质序列合理性
    
    📚 相关工具和应用 | Related Tools & Applications:
    
    🔧 上游工具:
    - 📊 基因组浏览器: IGV, UCSC Genome Browser
    - 🧬 基因组注释: NCBI, Ensembl, UCSC
    - 📄 格式转换: seqkit, seqtk, bioawk
    
    🔬 下游分析:
    - 🌳 序列比对: BLAST, MUSCLE, ClustalW
    - 🧪 PCR引物设计: Primer3, PrimerBLAST
    - 📊 系统发育分析: RAxML, IQ-TREE
    - 🔍 功能预测: InterProScan, SignalP
    
    📖 相关数据库:
    - 🌐 NCBI GenBank: 基因组序列数据
    - 🧬 Ensembl: 基因组注释数据
    - 📊 UCSC Genome Browser: 基因组区域信息
    - 🔍 UniProt: 蛋白质序列和功能
    """
    
    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    sequence_main = _lazy_import_sequence_main()
    
    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['parse_seq.py']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-s', sequence_file])
    args.extend(['-r', regions_file])
    args.extend(['-o', output_file])
    
    # 可选参数（只有在非默认值时才添加）⚙️ | Optional parameters (add only when non-default)
    if type != 'dna':
        args.extend(['--type', type])
    
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    if line_width != 80:
        args.extend(['--line-width', str(line_width)])
    
    if samtools_path != 'samtools':
        args.extend(['--samtools-path', samtools_path])
    
    # 处理互斥参数 | Handle mutually exclusive parameters
    if separate_output:
        args.append('--separate-output')
    elif not merge_output:
        # 如果merge_output被明确禁用但separate_output未启用
        pass
    
    if quiet:
        args.append('--quiet')
    elif not verbose:
        # 如果verbose被明确禁用但quiet未启用
        pass
    
    # 布尔选项 🚩 | Boolean options
    if no_headers:
        args.append('--no-headers')
    
    if reverse_complement:
        args.append('--reverse-complement')
    
    if translate:
        args.append('--translate')
    
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        sequence_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 序列提取被用户中断 | Sequence extraction interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 序列提取失败 | Sequence extraction failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv