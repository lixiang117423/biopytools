"""
VCF基因型提取命令 | VCF Genotype Extraction Command
"""

import click
import sys
from ...vcf_genotype_extractor.main import main as vcf_genotype_main


@click.command(short_help='VCF文件基因型数据提取和格式转换',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='VCF文件路径(支持.gz压缩格式) | VCF file path (supports .gz compressed format)')
@click.option('--output', '-o',
              default='vcf_genotype',
              type=str,
              help='输出文件前缀 | Output file prefix (default: vcf_genotype)')
@click.option('--samples', '-s',
              default='all',
              type=str,
              help='样本选择：all（所有样本）或逗号分隔的样本名称 | Sample selection: all (all samples) or comma-separated sample names (default: all)')
@click.option('--biallelic-only',
              is_flag=True,
              help='只保留双等位位点 | Keep only biallelic sites')
@click.option('--each', '-e',
              default='n',
              type=click.Choice(['yes', 'y', 'no', 'n']),
              help='按染色体拆分输出文件：yes/y（是）或no/n（否） | Split output files by chromosome: yes/y or no/n (default: n)')
@click.option('--output-type', '-t',
              default='txt',
              type=click.Choice(['txt', 'csv', 'excel']),
              help='输出文件格式 | Output file format (default: txt)')
@click.option('--output-dir',
              default='./',
              type=click.Path(),
              help='输出目录 | Output directory (default: ./)')
def vcf_genotype(input, output, samples, biallelic_only, each, output_type, output_dir):
    """
    VCF基因型提取工具
    
    从VCF文件中提取基因型数据，支持样本选择、格式转换和染色体分割输出。
    可处理大型VCF文件，支持压缩格式，提供多种输出选项。
    
    功能特点 | Features:
    - 高效VCF文件解析和基因型提取
    - 灵活的样本选择和过滤
    - 多种输出格式支持
    - 染色体级别数据分割
    - 压缩文件支持
    - 详细的处理统计
    
    输出格式 | Output Formats:
    - txt: 制表符分隔的文本文件
    - csv: 逗号分隔值文件
    - excel: Excel工作簿文件
    
    示例 | Examples:
    
    \b
    # 基本用法 - 提取所有样本基因型
    biopytools vcf-genotype -i variants.vcf -o genotypes
    
    \b
    # 提取特定样本
    biopytools vcf-genotype -i data.vcf.gz -o selected_samples \\
        -s "sample1,sample2,sample3"
    
    \b
    # 只保留双等位基因位点
    biopytools vcf-genotype -i variants.vcf -o biallelic_only \\
        --biallelic-only
    
    \b
    # 按染色体分别输出
    biopytools vcf-genotype -i genome.vcf.gz -o by_chromosome \\
        --each yes
    
    \b
    # Excel格式输出
    biopytools vcf-genotype -i data.vcf -o excel_output \\
        --output-type excel
    
    \b
    # 完整参数示例
    biopytools vcf-genotype -i large_cohort.vcf.gz \\
        -o cohort_genotypes -s "sample1,sample2,sample5" \\
        --biallelic-only --each yes --output-type csv \\
        --output-dir ./results/
    
    \b
    # 处理压缩文件并指定输出目录
    biopytools vcf-genotype -i compressed.vcf.gz \\
        -o extracted_data --output-dir /project/results/
    
    样本选择说明 | Sample Selection:
    - "all": 提取所有样本的基因型
    - "sample1,sample2": 提取指定样本（逗号分隔）
    - 工具会自动验证样本是否存在于VCF文件中
    
    输出文件命名 | Output File Naming:
    - 单文件模式: {prefix}.{format}
    - 染色体分割模式: {prefix}_{chromosome}.{format}
    - 汇总统计: {prefix}_summary.txt
    
    输出内容结构 | Output Content Structure:
    
    基因型矩阵格式:
    CHROM    POS    ID    REF    ALT    sample1    sample2    sample3
    chr1     1000   rs1   A      T      0/1        1/1        0/0
    chr1     2000   rs2   G      C      0/0        0/1        1/1
    
    基因型编码说明:
    - 0/0: 参考基因型纯合子
    - 0/1 或 1/0: 杂合子
    - 1/1: 变异基因型纯合子
    - ./.: 缺失数据
    
    处理流程 | Processing Pipeline:
    1. VCF文件格式检测和解析
    2. 样本信息验证和筛选
    3. 变异位点过滤（可选双等位基因）
    4. 基因型数据提取和格式化
    5. 输出文件生成（单文件或分染色体）
    6. 处理统计汇总报告
    
    应用场景 | Use Cases:
    - GWAS分析数据预处理
    - 群体遗传学矩阵构建
    - 系统发育分析数据准备
    - 基因型-表型关联研究
    - 育种数据管理
    
    性能优化 | Performance Features:
    - 内存高效的流式处理
    - 支持大型VCF文件
    - 压缩文件直接读取
    - 并行化数据处理
    - 智能缓存机制
    
    质量控制 | Quality Control:
    - 自动检测和处理缺失值
    - 基因型质量评估
    - 样本完整性检查
    - 详细的处理日志
    
    文件格式兼容性 | File Format Compatibility:
    - 标准VCF 4.0+格式
    - 压缩的.vcf.gz文件
    - bgzip压缩文件
    - 支持大型多样本VCF
    
    错误处理 | Error Handling:
    - 格式错误自动检测
    - 样本名称验证
    - 内存不足优雅处理
    - 详细错误报告
    
    输出统计信息 | Output Statistics:
    - 处理的变异位点数量
    - 提取的样本数量
    - 基因型完整性统计
    - 处理时间和性能指标
    
    注意事项 | Notes:
    - 大文件建议使用染色体分割输出
    - Excel格式适合小数据集可视化
    - CSV格式便于后续统计分析
    - 确保有足够磁盘空间存储结果
    """
    
    # 构建参数列表传递给原始main函数
    args = ['vcf_genotype.py']
    
    # 必需参数
    args.extend(['-i', input])
    
    # 可选参数（只在非默认值时添加）
    if output != 'vcf_genotype':
        args.extend(['-o', output])
    
    if samples != 'all':
        args.extend(['-s', samples])
    
    if biallelic_only:
        args.append('--biallelic-only')
    
    if each != 'n':
        args.extend(['-e', each])
    
    if output_type != 'txt':
        args.extend(['-t', output_type])
    
    if output_dir != './':
        args.extend(['--output-dir', output_dir])
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        vcf_genotype_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n提取流程被用户中断 | Extraction pipeline interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"提取失败 | Extraction failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv