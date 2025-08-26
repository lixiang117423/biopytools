"""
序列提取命令 | Sequence Extraction Command
"""

import click
import sys
from ...vcf_sequence_toolkit.main import main as sequence_extractor_main


@click.command(short_help='VCF和基因组序列变异提取工具',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--vcf', '-v',
              required=True,
              type=click.Path(exists=True),
              help='VCF文件路径 | VCF file path')
@click.option('--genome', '-g',
              required=True,
              type=click.Path(exists=True),
              help='基因组FASTA文件路径 | Genome FASTA file path')
@click.option('--chrom', '-c',
              required=True,
              type=str,
              help='染色体名称 | Chromosome name (e.g., chr1, 1)')
@click.option('--start', '-s',
              required=True,
              type=int,
              help='起始位置 (1-based) | Start position (1-based)')
@click.option('--end', '-e',
              required=True,
              type=int,
              help='结束位置 (1-based, inclusive) | End position (1-based, inclusive)')
@click.option('--output-dir', '-o',
              default='./sequence_output',
              type=click.Path(),
              help='输出目录 | Output directory (default: ./sequence_output)')
@click.option('--format',
              type=click.Choice(['tab', 'fasta', 'csv']),
              default='tab',
              help='输出格式 | Output format (default: tab)')
@click.option('--second-allele',
              is_flag=True,
              help='使用第二个等位基因而不是第一个 | Use second allele instead of first')
@click.option('--no-reference',
              is_flag=True,
              help='不包含参考序列 | Do not include reference sequence')
@click.option('--min-qual',
              type=int,
              help='最小质量值过滤 | Minimum quality filter')
@click.option('--samples',
              type=str,
              help='指定样品列表文件或逗号分隔的样品名称 | Sample list file or comma-separated sample names')
@click.option('--exclude-samples',
              type=str,
              help='排除样品列表文件或逗号分隔的样品名称 | Exclude sample list file or comma-separated sample names')
def vcf_sequence(vcf, genome, chrom, start, end, output_dir, format, second_allele,
                 no_reference, min_qual, samples, exclude_samples):
    """
    序列变异提取脚本 (模块化版本)
    
    从VCF文件和基因组文件中提取特定区间的序列变异信息，支持多种输出格式，
    适用于基因组学研究中的序列分析、变异检测验证和功能注释等应用场景。
    
    功能特点 | Features:
    - 🧬 精确的基因组区间序列提取
    - 📊 多种输出格式支持 (TAB, FASTA, CSV)
    - 🎯 灵活的样本选择和过滤
    - 🔍 变异质量控制和过滤
    - 📋 详细的统计报告生成
    - 💾 高效的大文件处理
    - 🚀 模块化可扩展架构
    
    分析流程 | Analysis Pipeline:
    1. 基因组和VCF文件验证
    2. 目标区间参考序列获取
    3. 区间内变异位点识别
    4. 样本基因型信息解析
    5. 个体化序列构建
    6. 多格式结果输出
    7. 统计报告生成
    
    应用场景 | Use Cases:
    - 🔬 候选基因序列分析
    - 🧪 功能验证实验序列准备
    - 📈 群体基因组比较分析
    - 🎯 特定区域变异模式研究
    - 💊 药物基因组学序列分析
    - 🌱 育种相关基因序列提取
    - 📊 进化和系统发育分析
    
    示例 | Examples:
    
    \b
    # 基本序列提取
    biopytools vcf-sequence -v variants.vcf -g genome.fa \\
        -c chr1 -s 1000000 -e 1001000 -o seq_output
    
    \b
    # 提取FASTA格式序列
    biopytools vcf-sequence -v variants.vcf.gz -g genome.fa \\
        -c 1 -s 1000000 -e 1001000 --format fasta -o fasta_output
    
    \b
    # 指定特定样本提取
    biopytools vcf-sequence -v population.vcf -g reference.fa \\
        -c chr2 -s 500000 -e 502000 --samples "Sample1,Sample2,Sample3"
    
    \b
    # 使用第二等位基因并应用质量过滤
    biopytools vcf-sequence -v variants.vcf -g genome.fa \\
        -c chr3 -s 2000000 -e 2005000 --second-allele --min-qual 30
    
    \b
    # 排除低质量样本并输出CSV格式
    biopytools vcf-sequence -v data.vcf -g genome.fa \\
        -c chrX -s 1000000 -e 1010000 --format csv \\
        --exclude-samples low_quality_samples.txt --min-qual 20
    
    \b
    # 大区间序列提取（不包含参考序列）
    biopytools vcf-sequence -v large_cohort.vcf.gz -g genome.fa \\
        -c chr22 -s 10000000 -e 15000000 --no-reference \\
        --format fasta -o large_region
    
    \b
    # 从样本列表文件提取特定群体序列
    biopytools vcf-sequence -v population_variants.vcf -g reference.fa \\
        -c chr7 -s 5000000 -e 5100000 --samples european_samples.txt \\
        --format tab --min-qual 25
    
    输入文件要求 | Input File Requirements:
    
    VCF文件格式:
    - 符合VCF 4.0+标准格式
    - 支持压缩格式 (.vcf.gz)
    - 包含完整的基因型信息
    - 建议包含质量分数 (QUAL字段)
    - 变异位点应覆盖目标提取区间
    
    基因组FASTA文件要求:
    - 标准FASTA格式
    - 染色体序列完整
    - 建议使用与VCF文件相同的参考基因组版本
    - 支持压缩格式 (.fa.gz, .fasta.gz)
    - 染色体命名应与VCF文件一致
    
    样本文件格式 (可选):
    - 每行一个样本名称
    - 样本名称应与VCF文件中的样本名匹配
    - 支持注释行 (以#开头)
    - 空行将被自动忽略
    
    参数详解 | Parameter Details:
    
    位置参数 | Position Parameters:
    --chrom (-c):
    - 指定目标染色体
    - 支持多种命名格式: "chr1", "1", "chrX"
    - 必须与VCF和基因组文件中的命名一致
    
    --start (-s) & --end (-e):
    - 使用1-based坐标系统
    - start和end位置都包含在提取区间内
    - 建议区间大小: 1bp - 10Mb
    - 超大区间可能影响处理速度和内存使用
    
    输出格式选项 | Output Format Options:
    
    --format tab (默认):
    - 制表符分隔的文本格式
    - 包含位置、参考序列、样本序列信息
    - 便于后续数据处理和分析
    - 适合大量样本的批量处理
    
    --format fasta:
    - 标准FASTA序列格式
    - 每个样本一个序列条目
    - 适合系统发育分析和序列比对
    - 序列ID包含样本名和区间信息
    
    --format csv:
    - 逗号分隔的CSV格式
    - 兼容Excel和其他电子表格软件
    - 适合统计分析和数据可视化
    - 包含完整的元数据信息
    
    等位基因选择 | Allele Selection:
    
    默认行为 (第一等位基因):
    - 对于杂合子基因型 (如0/1)，选择第一个等位基因
    - 适用于单倍体分析或随机等位基因选择
    
    --second-allele:
    - 对于杂合子基因型，选择第二个等位基因
    - 用于比较分析或特定等位基因研究
    - 与第一等位基因结果配合使用可进行单倍型分析
    
    质量控制选项 | Quality Control Options:
    
    --min-qual:
    - 设置变异位点的最小质量分数阈值
    - 过滤低质量的变异调用
    - 建议范围: 20-50
    - 提高序列构建的准确性
    
    --no-reference:
    - 排除参考序列的输出
    - 减少输出文件大小
    - 专注于样本间变异比较
    - 适用于已知参考序列的场景
    
    样本选择策略 | Sample Selection Strategy:
    
    --samples:
    - 正向选择: 仅处理指定的样本
    - 支持文件路径或逗号分隔的样本名
    - 文件格式: 每行一个样本名
    - 用于特定群体或感兴趣样本的分析
    
    --exclude-samples:
    - 反向选择: 排除指定的样本
    - 用于去除低质量样本或对照组
    - 与--samples参数互斥使用
    - 适用于大规模数据的质量控制
    
    输出文件说明 | Output Files Description:
    
    主要输出文件:
    - {chrom}_{start}_{end}_sequences.{format}: 提取的序列文件
    - {chrom}_{start}_{end}_statistics.txt: 提取统计信息
    - extraction_summary.txt: 提取过程总结
    - sequence_extraction.log: 详细的处理日志
    
    统计文件内容:
    - 提取区间基本信息
    - 参考序列长度和组成
    - 变异位点数量和类型分布
    - 样本序列质量统计
    - 等位基因频率分布
    
    日志文件记录:
    - 文件打开和验证过程
    - 序列提取的详细步骤
    - 变异位点处理信息
    - 错误和警告信息
    - 性能和时间统计
    
    技术细节 | Technical Details:
    
    内存管理:
    - 采用流式处理减少内存占用
    - 大区间自动分片处理
    - 智能缓存提高访问效率
    - 及时释放不再使用的资源
    
    性能优化:
    - 索引加速基因组访问
    - 并行处理多个样本序列
    - 高效的变异位点查找算法
    - 优化的文件I/O操作
    
    数据准确性:
    - 严格的坐标系统转换
    - 完整的基因型解析
    - 准确的序列重构算法
    - 全面的错误检测机制
    
    常见问题解决 | Troubleshooting:
    
    坐标系统问题:
    - 确保使用1-based坐标
    - 验证VCF和基因组文件的一致性
    - 检查染色体命名格式匹配
    
    内存不足:
    - 减小提取区间大小
    - 限制处理的样本数量
    - 使用--no-reference减少内存使用
    - 增加系统可用内存
    
    处理速度慢:
    - 检查磁盘I/O性能
    - 使用SSD存储提高速度
    - 考虑预先索引基因组文件
    - 优化VCF文件格式 (使用压缩)
    
    结果验证异常:
    - 检查参考基因组版本一致性
    - 验证VCF文件格式正确性
    - 确认目标区间包含足够的变异位点
    - 比对已知序列验证结果准确性
    
    最佳实践 | Best Practices:
    
    1. 数据准备:
       - 使用相同版本的参考基因组
       - 预先验证VCF文件完整性
       - 建立基因组文件索引加速访问
       - 确保染色体命名一致性
    
    2. 参数设置:
       - 根据研究目的选择合适的输出格式
       - 设定合理的质量阈值平衡质量与数量
       - 选择代表性样本避免冗余分析
       - 考虑下游分析需求设定输出选项
    
    3. 质量控制:
       - 验证提取区间的生物学意义
       - 检查变异密度是否符合预期
       - 比较不同样本间的序列差异
       - 验证关键位点的变异状态
    
    4. 结果应用:
       - 结合功能注释解释序列变异
       - 使用多种工具验证关键发现
       - 保留原始数据以备重复分析
       - 记录分析参数便于方法重现
    
    高级应用技巧 | Advanced Usage Tips:
    
    批量区间处理:
    - 编写脚本循环处理多个区间
    - 使用并行处理提高效率
    - 统一管理输出文件和目录
    
    单倍型分析:
    - 分别提取两个等位基因
    - 比较不同等位基因的序列差异
    - 结合连锁分析确定单倍型结构
    
    群体比较分析:
    - 按群体分别提取序列
    - 计算群体间序列差异
    - 识别群体特异性变异模式
    
    功能序列分析:
    - 提取编码序列进行翻译
    - 分析启动子区域的变异影响
    - 研究调控元件的序列保守性
    """
    
    # 构建参数列表传递给原始main函数
    args = ['vcf-sequence.py']
    
    # 必需参数
    args.extend(['-v', vcf])
    args.extend(['-g', genome])
    args.extend(['-c', chrom])
    args.extend(['-s', str(start)])
    args.extend(['-e', str(end)])
    
    # 可选参数（只在非默认值时添加）
    if output_dir != './sequence_output':
        args.extend(['-o', output_dir])
    
    if format != 'tab':
        args.extend(['--format', format])
    
    if min_qual is not None:
        args.extend(['--min-qual', str(min_qual)])
    
    if samples:
        args.extend(['--samples', samples])
    
    if exclude_samples:
        args.extend(['--exclude-samples', exclude_samples])
    
    # 布尔选项
    if second_allele:
        args.append('--second-allele')
    
    if no_reference:
        args.append('--no-reference')
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        result = sequence_extractor_main()
        sys.exit(result if result is not None else 0)
    except SystemExit as e:
        # 处理程序正常退出
        sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n序列提取被用户中断 | Sequence extraction interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"序列提取失败 | Sequence extraction failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv