"""
VCF转换工具命令 | VCF Converter Command
"""

import click
from ...vcf2phylip.main import VCFConverter


@click.command(short_help='VCF格式转换工具',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='输入VCF文件路径 (支持gzip压缩) | Input VCF file path (gzip supported)')
@click.option('--output', '-o',
              default='./converted_output',
              type=click.Path(),
              help='输出目录 | Output directory (default: ./converted_output)')
@click.option('--output-prefix',
              type=str,
              help='输出文件名前缀 | Output filename prefix')
@click.option('--min-samples-locus', '-m',
              default=4,
              type=int,
              help='位点最少样本数 | Minimum samples required per locus (default: 4)')
@click.option('--outgroup', '-g',
              default="",
              type=str,
              help='外群样本名称 | Outgroup sample name (default: "")')
@click.option('--phylip-disable', '-p',
              is_flag=True,
              help='禁用PHYLIP输出 | Disable PHYLIP output')
@click.option('--fasta', '-f',
              is_flag=True,
              help='启用FASTA输出 | Enable FASTA output')
@click.option('--nexus', '-n',
              is_flag=True,
              help='启用NEXUS输出 | Enable NEXUS output')
@click.option('--nexus-binary', '-b',
              is_flag=True,
              help='启用二进制NEXUS输出 (仅适用于二倍体) | Enable binary NEXUS output (diploid only)')
@click.option('--resolve-IUPAC', '-r',
              is_flag=True,
              help='随机解析杂合子基因型以避免IUPAC模糊代码 | Randomly resolve heterozygous genotypes')
@click.option('--write-used-sites', '-w',
              is_flag=True,
              help='保存通过筛选的位点坐标列表 | Save list of coordinates that passed filters')
@click.option('--threads', '-t',
              default=88,
              type=int,
              help='线程数 | Number of threads (default: 88)')
@click.version_option(version='v2.9.1', message='%(prog)s %(version)s - 模块化版本 | Modular Version')
def vcf2phylip(input, output, output_prefix, min_samples_locus, outgroup,
                  phylip_disable, fasta, nexus, nexus_binary, resolve_iupac,
                  write_used_sites, threads):
    """
    VCF转换工具 (模块化版本)
    
    将VCF格式的SNP数据转换为多种系统发育分析格式，支持PHYLIP、FASTA、NEXUS
    和二进制NEXUS格式，适用于群体基因组学和系统发育研究。
    
    功能特点 | Features:
    - 多种输出格式支持 (PHYLIP/FASTA/NEXUS)
    - 高效的多线程处理
    - 灵活的样本和位点过滤
    - IUPAC模糊代码处理
    - 外群样本指定
    - 二进制NEXUS格式支持
    - 大文件优化处理
    
    转换流程 | Conversion Pipeline:
    1. VCF文件格式验证和解析
    2. 样本名称提取和验证
    3. SNP位点质量过滤
    4. 基因型数据处理
    5. 多格式矩阵生成
    6. 外群处理和排序
    7. 结果文件输出
    
    应用场景 | Use Cases:
    - 系统发育树构建预处理
    - 群体结构分析数据准备
    - 比较基因组学研究
    - 进化分析输入格式转换
    - 大规模SNP数据处理
    - 多物种基因组比较
    
    示例 | Examples:
    
    \b
    # 基本VCF转换 (默认PHYLIP格式)
    biopytools vcf-converter -i variants.vcf -o converted_results
    
    \b
    # 转换为多种格式
    biopytools vcf-converter -i data.vcf.gz -o results --fasta --nexus
    
    \b
    # 指定外群和最小样本数
    biopytools vcf-converter -i variants.vcf -o results \\
        -m 10 -g outgroup_sample
    
    \b
    # 二进制NEXUS格式转换
    biopytools vcf-converter -i data.vcf -o results \\
        --nexus-binary --resolve-IUPAC
    
    \b
    # 高性能批量转换
    biopytools vcf-converter -i large_dataset.vcf.gz \\
        -o batch_results --threads 64 --fasta --nexus
    
    \b
    # 完整转换分析
    biopytools vcf-converter -i population.vcf -o comprehensive \\
        --fasta --nexus --nexus-binary -m 8 -g outgroup \\
        --resolve-IUPAC --write-used-sites --threads 32
    
    输入文件要求 | Input File Requirements:
    
    VCF文件格式:
    - 符合VCF 4.0+标准
    - 支持压缩格式 (.vcf.gz)
    - 包含完整的基因型信息 (GT字段)
    - 建议包含质量信息 (QUAL字段)
    - 样本名称不能重复
    
    数据质量要求:
    - SNP位点应有足够的样本覆盖
    - 基因型调用质量较高
    - 避免过多的缺失数据
    - 建议预先进行基本质控
    
    输出格式说明 | Output Format Details:
    
    PHYLIP格式 (默认):
    - 适用于RAxML, IQ-TREE等软件
    - 严格PHYLIP格式或relaxed格式
    - 包含样本数和位点数信息
    - 序列名称长度限制处理
    
    FASTA格式:
    - 标准FASTA序列格式
    - 适用于多种系统发育软件
    - 样本名作为序列标识符
    - 支持长序列名称
    
    NEXUS格式:
    - 适用于MrBayes, BEAST等软件
    - 包含完整的数据块定义
    - 支持字符状态定义
    - 可包含外群信息
    
    二进制NEXUS格式:
    - 仅适用于二倍体数据
    - 将SNP编码为0/1二进制
    - 适用于特定的贝叶斯分析
    - 显著减少文件大小
    
    参数详解 | Parameter Details:
    
    --min-samples-locus (-m):
    - 设置每个SNP位点的最小样本数
    - 过滤覆盖度不足的位点
    - 建议范围: 3-样本总数的80%
    - 平衡数据质量与数据量
    
    --outgroup (-g):
    - 指定外群样本名称
    - 外群将排列在输出矩阵首位
    - 样本名必须完全匹配VCF中的名称
    - 用于后续系统发育分析的根定位
    
    --resolve-IUPAC (-r):
    - 随机解析杂合子位点
    - 避免IUPAC模糊碱基代码
    - 适用于不支持模糊代码的软件
    - 每次运行结果可能不同
    
    --write-used-sites (-w):
    - 输出通过筛选的位点坐标
    - 便于后续分析的位点追踪
    - 生成sites_used.txt文件
    - 包含染色体和位置信息
    
    性能优化建议 | Performance Optimization:
    
    线程配置:
    - 小文件 (<10MB): 4-8线程
    - 中等文件 (10MB-1GB): 16-32线程
    - 大文件 (>1GB): 32-88线程
    
    内存管理:
    - 大VCF文件使用分块处理
    - 临时文件自动管理
    - 内存使用随样本数和位点数增长
    
    输出格式选择:
    - PHYLIP: 最通用，文件较小
    - FASTA: 兼容性好，易于处理
    - NEXUS: 功能丰富，文件较大
    - 二进制NEXUS: 特殊用途，最小文件
    
    质量控制策略 | Quality Control Strategy:
    
    输入验证:
    - 自动检测VCF格式完整性
    - 验证样本名称唯一性
    - 检查基因型字段完整性
    
    数据过滤:
    - 基于样本覆盖度过滤位点
    - 自动调整参数以适应数据
    - 保留高质量SNP位点
    
    输出验证:
    - 检查输出文件格式正确性
    - 验证序列长度一致性
    - 统计转换结果信息
    
    常见问题解决 | Troubleshooting:
    
    内存不足:
    - 减少线程数
    - 分批处理大文件
    - 关闭不需要的输出格式
    
    格式错误:
    - 检查VCF文件完整性
    - 验证样本名称格式
    - 确认基因型字段存在
    
    转换失败:
    - 检查输入文件权限
    - 验证输出目录可写
    - 查看详细错误日志
    
    性能慢:
    - 适当增加线程数
    - 使用SSD存储
    - 预处理VCF文件质量
    
    外群问题:
    - 确保外群样本名完全匹配
    - 检查外群样本在VCF中存在
    - 验证外群数据质量
    
    高级应用技巧 | Advanced Usage Tips:
    
    批量处理:
    - 编写脚本处理多个VCF文件
    - 使用一致的参数设置
    - 统一管理输出文件命名
    
    数据预处理:
    - 使用VCFtools进行预过滤
    - 去除低质量样本和位点
    - 标准化VCF格式
    
    下游分析集成:
    - 直接对接系统发育分析软件
    - 保持参数设置的一致性
    - 记录转换参数用于重现
    
    特殊数据处理:
    - 古DNA数据: 适当降低样本数要求
    - 多倍体数据: 避免使用二进制格式
    - 大规模数据: 考虑分染色体处理
    """
    
    try:
        # 直接创建转换器，避免参数重构的性能损失
        converter = VCFConverter(
            input_file=input,
            output_dir=output,
            output_prefix=output_prefix,
            min_samples_locus=min_samples_locus,
            outgroup=outgroup,
            phylip_disable=phylip_disable,
            fasta=fasta,
            nexus=nexus,
            nexus_binary=nexus_binary,
            resolve_IUPAC=resolve_iupac,
            write_used_sites=write_used_sites,
            threads=threads
        )
        
        # 直接运行转换，无额外包装
        converter.run_conversion()
        
    except KeyboardInterrupt:
        click.echo("\nVCF转换被用户中断 | VCF conversion interrupted by user", err=True)
        raise click.Abort()
    except Exception as e:
        click.echo(f"VCF转换失败 | VCF conversion failed: {e}", err=True)
        raise click.Abort()