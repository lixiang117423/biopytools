"""
GenBank序列提取命令 | GenBank Sequence Extraction Command
高级优化版本：解决--help响应速度问题
"""

import click
import sys
import os


def _lazy_import_genbank_main():
    """懒加载GenBank main函数 | Lazy load GenBank main function"""
    try:
        from ...genebank2fasta.main import main as genbank_main
        return genbank_main
    except ImportError as e:
        click.echo(f"❌ 导入错误 | Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否是帮助请求 | Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_directory_exists(dir_path):
    """验证目录是否存在（仅在非帮助模式下）| Validate directory existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(dir_path):
        raise click.BadParameter(f"❌ 目录不存在 | Directory does not exist: {dir_path}")
    if not _is_help_request() and not os.path.isdir(dir_path):
        raise click.BadParameter(f"❌ 路径不是目录 | Path is not a directory: {dir_path}")
    return dir_path


@click.command(short_help='🧬 GenBank序列提取工具：从GenBank文件批量提取CDS和蛋白质序列',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_directory_exists(value) if value else None,
              help='📂 输入GenBank文件目录 | Input GenBank files directory')
@click.option('--output', '-o',
              default='./genbank_output',
              type=click.Path(),
              help='📁 输出目录 | Output directory (default: ./genbank_output)')
@click.option('--threads', '-t',
              type=int,
              default=88,
              help='⚡ 并行线程数 | Number of parallel threads (default: 88)')
@click.option('--min-length',
              type=int,
              default=10,
              help='🔬 最小蛋白质长度 (氨基酸) | Minimum protein length (amino acids) (default: 10)')
@click.option('--phylo', '--create-phylogenetic-matrix',
              is_flag=True,
              help='🌲 创建系统发育分析矩阵 | Create phylogenetic analysis matrix')
@click.option('--no-sample-sep',
              is_flag=True,
              help='🚫 不按样品分离输出 | Do not separate output by sample')
@click.option('--no-gene-sep',
              is_flag=True,
              help='🚫 不按基因分离输出 | Do not separate output by gene')
@click.option('--keep-unknown',
              is_flag=True,
              help='📝 保留unknown基因 | Keep unknown genes')
def genebank2fasta(input, output, threads, min_length, phylo, no_sample_sep, 
                     no_gene_sep, keep_unknown):
    """
    🧬 GenBank序列提取工具 | GenBank Sequence Extraction Tool
    
    一个专业的GenBank格式文件批量处理工具，用于从GenBank (.gb/.gbk/.genbank)文件中
    自动提取编码序列(CDS)和对应的蛋白质序列，支持大规模并行处理和多种输出格式，
    适用于比较基因组学、系统发育分析和功能基因组学研究。
    
    ✨ 功能特点 | Features:
    - 🚀 高效的并行批量处理能力
    - 📊 智能的GenBank文件格式解析
    - 🧬 同时提取CDS核酸和蛋白质序列
    - 🏷️ 灵活的序列命名和组织方式
    - 🌲 可选的系统发育分析矩阵生成
    - 📈 详细的提取统计和质量报告
    - 🔍 可定制的序列长度过滤
    
    🔄 处理流程 | Processing Pipeline:
    1. 📂 扫描输入目录中的GenBank文件
    2. 🏗️ 创建结构化的输出目录体系
    3. ⚡ 多线程并行解析GenBank文件
    4. 🧬 提取CDS特征和翻译蛋白质序列
    5. 📊 按基因或样本分类组织序列
    6. 💾 生成FASTA格式的输出文件
    7. 📈 生成综合统计报告
    8. 🌲 构建系统发育分析矩阵(可选)
    
    🎯 应用场景 | Use Cases:
    - 🔬 比较基因组学序列数据准备
    - 🌳 系统发育分析序列收集
    - 🧪 功能基因组学研究数据处理
    - 📊 基因家族分析序列整理
    - 🧬 进化分析数据预处理
    - 💾 基因组注释数据提取
    - 🔍 特定基因的批量序列获取
    
    💡 示例 | Examples:
    
    \b
    # 🎯 基本序列提取
    biopytools genbank-extractor -i /path/to/genbank/files -o ./output
    
    \b
    # ⚡ 高性能并行处理
    biopytools genbank-extractor -i ./gb_files -o ./results -t 64
    
    \b
    # 🌲 包含系统发育矩阵生成
    biopytools genbank-extractor -i ./genbank -o ./extract_results \\
        -t 32 --phylo
    
    \b
    # 🔬 自定义最小长度和输出组织
    biopytools genbank-extractor -i /data/gb_files -o /results \\
        --min-length 50 --threads 88 --no-sample-sep
    
    \b
    # 📝 保留所有基因包括未知功能
    biopytools genbank-extractor -i ./bacterial_genomes -o ./all_genes \\
        --keep-unknown --no-gene-sep
    
    \b
    # 🧬 系统发育研究专用设置
    biopytools genbank-extractor -i ./genomes_collection -o ./phylo_data \\
        --phylo --min-length 100 -t 96 --no-sample-sep
    
    📁 输入文件要求 | Input File Requirements:
    
    🗂️ GenBank文件格式:
    - 📄 标准GenBank格式文件 (.gb, .gbk, .genbank)
    - 🗜️ 支持压缩文件 (.gb.gz, .gbk.gz)
    - 🧬 包含完整的CDS特征注释
    - 🏷️ 包含基因名称和产品描述
    - ✅ 符合NCBI GenBank格式规范
    
    📂 目录结构要求:
    - 📁 所有GenBank文件放在指定输入目录中
    - 🔍 支持递归搜索子目录中的文件
    - 📝 文件命名建议包含样本标识信息
    - 🗃️ 可包含混合的不同物种基因组文件
    
    📄 GenBank格式示例:
    LOCUS       NC_000001    50000 bp    DNA     linear   CON 01-JAN-2023
    DEFINITION  Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly
    ...
    FEATURES             Location/Qualifiers
         CDS             1000..2500
                         /gene="example_gene"
                         /product="hypothetical protein"
                         /translation="MKLAVF..."
    
    ⚙️ 参数详细说明 | Parameter Details:
    
    📋 必需参数:
    --input (-i): 📂 GenBank文件输入目录路径
        - 包含.gb/.gbk/.genbank文件的目录
        - 支持递归搜索子目录
        - 自动识别压缩文件
    
    📁 输出控制:
    --output (-o): 📁 结果输出目录路径
        - 自动创建必要的子目录结构
        - 默认在当前目录下创建genbank_output文件夹
    
    ⚡ 性能参数:
    --threads (-t): 🧵 并行处理线程数
        - 建议设置为CPU核心数
        - 影响文件处理和序列提取速度
        - 大文件数量时显著提升性能
    
    🔬 质量控制:
    --min-length: 📏 最小蛋白质序列长度
        - 单位为氨基酸残基数量
        - 过滤过短的可能错误的CDS
        - 有助于提高序列质量
    
    📊 输出组织选项:
    --no-sample-sep: 🚫 不按样本分离输出
        - 将所有样本的序列合并到单一文件
        - 适用于基因家族分析
        - 减少输出文件数量
    
    --no-gene-sep: 🚫 不按基因分离输出
        - 将所有基因序列合并输出
        - 适用于全基因组序列分析
        - 生成样本级别的序列文件
    
    --keep-unknown: 📝 保留未知功能基因
        - 包含product为"hypothetical protein"的基因
        - 保留功能注释缺失的CDS
        - 适用于全基因组分析研究
    
    🌲 高级功能:
    --phylo: 🌳 创建系统发育分析矩阵
        - 生成适合系统发育分析的序列矩阵
        - 自动处理序列比对格式
        - 输出多种系统发育软件兼容格式
    
    📁 输出文件结构 | Output File Structure:
    
    📂 标准输出目录结构:
    output_directory/
    ├── 📁 cds/                    # CDS核酸序列
    │   ├── 📄 by_gene/           # 按基因分类的CDS文件
    │   └── 📄 by_sample/         # 按样本分类的CDS文件
    ├── 📁 proteins/               # 蛋白质序列
    │   ├── 📄 by_gene/           # 按基因分类的蛋白质文件
    │   └── 📄 by_sample/         # 按样本分类的蛋白质文件
    ├── 📊 statistics/             # 统计报告
    │   ├── 📈 extraction_summary.txt
    │   ├── 📊 gene_statistics.txt
    │   └── 📋 sample_statistics.txt
    └── 🌲 phylogenetic/           # 系统发育分析矩阵(可选)
        ├── 📄 concatenated_cds.fasta
        ├── 📄 concatenated_proteins.fasta
        └── 📋 gene_positions.txt
    
    📄 序列文件命名规则:
    按基因分类:
    - 🧬 gene_name.cds.fasta: 特定基因的CDS序列
    - 🧪 gene_name.pep.fasta: 特定基因的蛋白质序列
    
    按样本分类:
    - 📊 sample_name.cds.fasta: 样本的所有CDS序列
    - 📊 sample_name.pep.fasta: 样本的所有蛋白质序列
    
    📈 统计报告内容:
    - 📊 每个文件的处理结果统计
    - 🧬 基因提取数量和成功率
    - 📏 序列长度分布统计
    - ❌ 处理失败文件的详细信息
    - 🔍 质量过滤统计信息
    
    ⚡ 性能和系统要求 | Performance & System Requirements:
    
    💻 系统建议:
    - 💾 RAM: 至少4GB，大数据集建议16GB+
    - 🗄️ 存储: 输出目录需要输入数据2-3倍空间
    - 🧵 CPU: 多核处理器显著提升并行处理效率
    - 🐍 Python: 3.7+版本
    
    📊 处理能力估算:
    - 📄 小数据集(<100个文件): 几分钟内完成
    - 📋 中等数据集(100-1000个文件): 10-30分钟
    - 📁 大数据集(1000+个文件): 30分钟-2小时
    - 🧵 多线程可线性提升处理速度
    
    💾 内存使用估算:
    - 🔵 基础内存: ~500MB
    - 📈 文件相关: 每个大文件额外~50-100MB
    - 🧵 线程相关: 每线程额外~20MB
    
    🛠️ 故障排除 | Troubleshooting:
    
    ⚠️ 常见问题:
    1. ❌ "No GenBank files found": 检查输入目录路径和文件扩展名
    2. 💾 "内存不足": 减少线程数或处理更小的文件批次
    3. 📄 "文件格式错误": 验证GenBank文件格式完整性
    4. 🔒 "权限拒绝": 检查输出目录写入权限
    5. 🧬 "CDS特征缺失": 确认GenBank文件包含完整注释
    
    📋 数据质量检查:
    - 🔍 验证GenBank文件完整性和格式
    - 📏 检查CDS序列长度合理性
    - 🧪 确认蛋白质翻译正确性
    - 🏷️ 验证基因名称和产品注释
    
    💡 优化建议:
    - 📊 根据文件数量合理设置线程数
    - 💾 使用高速存储提升I/O性能
    - 🔍 预先检查文件质量避免处理失败
    - 📈 监控系统资源使用情况
    
    🏆 最佳实践 | Best Practices:
    
    1️⃣ 📊 数据准备:
       - ✅ 使用标准的GenBank格式文件
       - 🔍 预先验证文件完整性和可读性
       - 📝 使用有意义的文件命名规范
       - 🗂️ 按项目或研究目标组织文件
    
    2️⃣ ⚙️ 参数选择:
       - 🧵 根据系统配置设置合适线程数
       - 📏 基于研究需求调整最小长度阈值
       - 📊 考虑下游分析选择输出组织方式
       - 🌲 大规模系统发育研究启用矩阵功能
    
    3️⃣ 🔍 质量控制:
       - 📈 检查统计报告中的成功率
       - 🎲 随机验证提取序列的正确性
       - 📊 比较预期与实际的基因数量
       - 🧬 验证CDS和蛋白质序列一致性
    
    4️⃣ 📦 结果管理:
       - 💾 及时备份重要的提取结果
       - 📋 保留详细的处理日志记录
       - 🗂️ 建立清晰的文件组织结构
       - 📊 记录处理参数便于重现分析
    
    🔍 结果验证建议 | Result Validation:
    
    🤖 自动验证:
    - 📊 统计报告数据完整性检查
    - 🧬 序列格式和长度分布验证
    - 🏷️ 基因名称和注释信息核实
    - 📈 提取成功率和错误率统计
    
    👥 人工验证:
    - 🎲 随机抽查序列提取准确性
    - 🔍 验证关键基因的序列完整性
    - 📚 比较已知基因与提取结果
    - 🌐 与公共数据库序列进行比对
    
    📚 相关工具和数据库 | Related Tools & Databases:
    
    🔧 上游工具:
    - 📊 NCBI Genome下载工具
    - 🧬 基因组注释软件(Prokka, RAST等)
    - 📄 格式转换工具(seqret, convertio等)
    
    🔬 下游分析工具:
    - 🌳 系统发育分析: RAxML, IQ-TREE, MrBayes
    - 🔍 序列比对: MUSCLE, ClustalW, MAFFT
    - 📊 比较基因组学: OrthoFinder, GET_HOMOLOGUES
    - 🧪 功能注释: InterProScan, BLAST, HMMer
    
    📖 相关数据库:
    - 🌐 NCBI GenBank: 主要的基因组数据源
    - 🧬 RefSeq: 参考序列数据库
    - 🔍 UniProt: 蛋白质序列和功能数据库
    - 📊 KEGG: 代谢通路和功能分类数据库
    """
    
    # 🚀 懒加载：只有在实际调用时才导入模块 | Lazy loading: import only when actually called
    genbank_main = _lazy_import_genbank_main()
    
    # 构建参数列表传递给原始main函数 🔄 | Build argument list for original main function
    args = ['genebank2fasta.py']
    
    # 必需参数 📋 | Required parameters
    args.extend(['-i', input])
    
    # 可选参数（只有在非默认值时才添加）⚙️ | Optional parameters (add only when non-default)
    if output != './genbank_output':
        args.extend(['-o', output])
    
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    if min_length != 10:
        args.extend(['--min-length', str(min_length)])
    
    # 布尔选项 🚩 | Boolean options
    if phylo:
        args.append('--phylo')
    
    if no_sample_sep:
        args.append('--no-sample-sep')
    
    if no_gene_sep:
        args.append('--no-gene-sep')
    
    if keep_unknown:
        args.append('--keep-unknown')
    
    # 保存并恢复sys.argv 💾 | Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数 🚀 | Call original main function
        genbank_main()
    except SystemExit as e:
        # 处理程序正常退出 ✅ | Handle normal program exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n🛑 GenBank序列提取被用户中断 | GenBank sequence extraction interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 GenBank序列提取失败 | GenBank sequence extraction failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv