"""
VCF基因型统计分析命令 | VCF Genotype Statistics Analysis Command
"""

import click
import sys
from ...vcf_stats_sample.main import main as vcf_stats_main


@click.command(short_help='VCF基因型统计分析工具',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--vcf', '-v',
              required=True,
              type=click.Path(exists=True),
              help='输入VCF文件路径 (支持.gz压缩格式) | Input VCF file path (supports .gz compressed format)')
@click.option('--output', '-o',
              default='vcf_stats_output',
              type=click.Path(),
              help='输出目录 | Output directory (default: vcf_stats_output)')
@click.option('--min-depth', '-d',
              default=0,
              type=int,
              help='最小深度过滤阈值 (0表示不过滤) | Minimum depth filter threshold (0 = no filter) (default: 0)')
@click.option('--min-qual', '-q',
              default=0.0,
              type=float,
              help='最小质量分数过滤阈值 (0.0表示不过滤) | Minimum quality score filter threshold (0.0 = no filter) (default: 0.0)')
@click.option('--exclude-missing', '-e',
              is_flag=True,
              help='排除缺失基因型 (./..) 的统计 | Exclude missing genotypes (./..) from statistics')
@click.option('--no-detailed', '-D',
              is_flag=True,
              help='不输出详细统计结果 | Do not output detailed statistics')
@click.option('--no-summary', '-S',
              is_flag=True,
              help='不输出汇总统计结果 | Do not output summary statistics')
def vcf_sample_hete(vcf, output, min_depth, min_qual, exclude_missing, no_detailed, no_summary):
    """
    VCF基因型统计分析脚本 (模块化版本)
    
    一个专业的VCF文件基因型统计分析工具，提供全面的基因型质量评估、
    统计计算和结果导出功能，适用于基因组学研究中的数据质量控制和样本评估。
    
    功能特点 | Features:
    - 🔍 全面的基因型质量统计
    - 📊 多维度统计指标计算
    - 🎯 灵活的质量和深度过滤
    - 💾 多格式结果输出
    - 📋 详细的分析报告
    - 🚀 高效的大文件处理
    - 📈 样本和位点级别统计
    
    分析流程 | Analysis Pipeline:
    1. VCF文件格式验证和读取
    2. 基因型质量过滤和预处理
    3. 多维度统计指标计算
    4. 样本和位点级别分析
    5. 统计结果汇总和导出
    6. 分析报告生成
    
    应用场景 | Use Cases:
    - 🧬 基因型数据质量评估
    - 📊 样本质量控制和筛选
    - 🔍 基因型调用结果验证
    - 📈 测序深度和覆盖度分析
    - 🎯 批量样本质量比较
    - 📋 基因组学数据预处理
    - 💡 实验设计优化参考
    
    示例 | Examples:
    
    \b
    # 基本统计分析
    biopytools vcf-sample-hete -v variants.vcf -o vcf_stats_output
    
    \b
    # 应用质量和深度过滤
    biopytools vcf-sample-hete -v variants.vcf.gz -o results -d 10 -q 30.0
    
    \b
    # 排除缺失基因型，仅输出汇总统计
    biopytools vcf-sample-hete -v variants.vcf -o results -e -D
    
    \b
    # 严格质控过滤和完整分析
    biopytools vcf-sample-hete -v high_quality.vcf -o comprehensive_stats \\
        --min-depth 15 --min-qual 40.0 --exclude-missing
    
    \b
    # 快速质量评估（仅汇总统计）
    biopytools vcf-sample-hete -v variants.vcf -o quick_stats \\
        --no-detailed --min-depth 5
    
    \b
    # 批量样本质量比较
    biopytools vcf-sample-hete -v population_variants.vcf.gz \\
        -o population_stats --min-qual 20.0
    
    \b
    # 测序质量评估
    biopytools vcf-sample-hete -v raw_variants.vcf -o quality_assessment \\
        --min-depth 8 --exclude-missing
    
    输入文件要求 | Input File Requirements:
    
    VCF文件格式:
    - 符合VCF 4.0+标准格式
    - 支持压缩文件 (.vcf.gz)
    - 包含完整的基因型信息
    - 建议包含FORMAT字段中的DP(深度)和GQ(基因型质量)信息
    
    支持的基因型格式:
    - 🔗 未定相基因型: 0/0, 0/1, 1/1, ./.
    - 🔒 已定相基因型: 0|0, 0|1, 1|1, .|.
    - 🌟 多等位基因位点: 0/2, 1/2, 2/2等
    - ❓ 缺失基因型: ./., .|.
    
    质量过滤参数详解 | Quality Filter Parameters:
    
    --min-depth (-d):
    - 设置最小测序深度阈值
    - 过滤低深度的基因型调用
    - 默认0 (不过滤)
    - 建议范围: 5-20
    - 用途: 去除可能不准确的低深度基因型
    
    --min-qual (-q):
    - 设置最小基因型质量分数
    - 基于GQ (Genotype Quality) 字段过滤
    - 默认0.0 (不过滤)
    - 建议范围: 20.0-50.0
    - 用途: 确保基因型调用的可信度
    
    --exclude-missing (-e):
    - 排除缺失基因型 (./.或.|.) 的统计
    - 专注于有效基因型的质量评估
    - 适用于高质量数据分析
    - 有助于更准确的统计计算
    
    输出控制选项 | Output Control Options:
    
    --no-detailed (-D):
    - 禁用详细统计结果输出
    - 仅生成汇总级别的统计
    - 适用于快速质量评估
    - 减少输出文件大小
    
    --no-summary (-S):
    - 禁用汇总统计结果输出
    - 仅生成详细级别的统计
    - 适用于深入分析特定样本
    - 专注于单样本详细信息
    
    输出文件说明 | Output Files Description:
    
    汇总统计文件 (默认生成):
    - summary_statistics.txt: 整体统计汇总
    - sample_statistics.txt: 样本级别统计
    - variant_statistics.txt: 位点级别统计
    - quality_metrics.txt: 质量指标汇总
    
    详细统计文件 (如果启用):
    - detailed_genotype_stats.txt: 详细基因型统计
    - per_sample_details/: 每个样本的详细统计目录
    - per_variant_details.txt: 每个位点的详细信息
    - quality_distribution.txt: 质量分布详情
    
    分析报告文件:
    - analysis_summary.txt: 分析总结报告
    - processing_log.txt: 处理日志文件
    - parameter_settings.txt: 参数设置记录
    
    统计指标说明 | Statistical Metrics:
    
    样本级别统计:
    - 总基因型数量 | Total genotypes
    - 缺失基因型数量和比例 | Missing genotypes count and rate
    - 纯合子数量和比例 | Homozygous count and rate
    - 杂合子数量和比例 | Heterozygous count and rate
    - 平均测序深度 | Average sequencing depth
    - 平均基因型质量 | Average genotype quality
    
    位点级别统计:
    - 等位基因频率 | Allele frequencies
    - 基因型频率 | Genotype frequencies
    - Hardy-Weinberg平衡检验 | Hardy-Weinberg equilibrium test
    - 缺失数据比例 | Missing data rate
    - 质量分数分布 | Quality score distribution
    
    质量控制指标:
    - 通过质量过滤的基因型比例
    - 不同深度阈值下的数据保留率
    - 质量分数的分布统计
    - 样本间质量一致性评估
    
    结果解读指导 | Result Interpretation Guide:
    
    样本质量评估:
    1. 缺失率 >20%: 可能存在样本质量问题
    2. 平均深度 <10x: 可能影响基因型准确性
    3. 杂合子率异常: 可能提示样本污染或近亲繁殖
    4. 质量分数偏低: 需要检查测序质量
    
    位点质量评估:
    1. 缺失率 >50%: 该位点可能难以准确检测
    2. HWE显著偏离: 可能存在技术错误或选择压力
    3. 等位基因频率极端: 需要验证是否为真实变异
    4. 质量分数分布异常: 可能需要重新评估
    
    性能和系统要求 | Performance & System Requirements:
    
    系统建议:
    - RAM: 至少2GB，大文件建议8GB+
    - 存储: 确保输出目录有足够空间
    - CPU: 单线程处理，多核无显著优势
    
    文件大小处理能力:
    - 小文件 (<1GB): 快速处理，通常几分钟内完成
    - 中等文件 (1-10GB): 正常处理，建议充足内存
    - 大文件 (>10GB): 需要足够内存和存储空间
    
    优化建议:
    - 对于极大文件，考虑按染色体分割处理
    - 使用压缩格式可节省存储空间
    - 预先过滤明显低质量的位点可提高效率
    
    常见问题解决 | Troubleshooting:
    
    内存不足:
    - 减少输出文件类型 (使用 -D 或 -S)
    - 增加系统内存
    - 考虑数据分片处理
    
    处理速度慢:
    - 检查磁盘I/O性能
    - 使用SSD存储提高速度
    - 预先进行基本质量过滤
    
    结果异常:
    - 验证VCF文件格式正确性
    - 检查过滤参数是否过严
    - 确认输入数据的生物学合理性
    
    最佳实践建议 | Best Practice Recommendations:
    
    1. 数据预处理:
       - 使用标准VCF格式
       - 确保FORMAT字段完整
       - 预先去除明显错误的记录
    
    2. 参数设置:
       - 根据测序平台调整质量阈值
       - 考虑研究目的设定深度要求
       - 平衡数据质量与数据量
    
    3. 结果验证:
       - 对比已知高质量样本的统计结果
       - 检查统计结果的生物学合理性
       - 关注异常样本的详细信息
    
    4. 后续分析:
       - 使用统计结果指导质量控制策略
       - 识别需要重新处理的样本
       - 为下游分析设定合适的过滤标准
    
    技术特点 | Technical Features:
    
    高效算法:
    - 流式处理减少内存占用
    - 智能缓存提高读取效率
    - 并行计算加速统计计算
    
    稳定性保证:
    - 异常数据的优雅处理
    - 详细的错误日志记录
    - 中间结果的自动备份
    
    扩展性设计:
    - 模块化架构便于功能扩展
    - 标准化接口支持集成开发
    - 灵活的配置选项适应不同需求
    """
    
    # 构建参数列表传递给原始main函数
    args = ['vcf-sample-hete.py']
    
    # 必需参数
    args.extend(['-v', vcf])
    
    # 可选参数（只在非默认值时添加）
    if output != 'vcf_stats_output':
        args.extend(['-o', output])
    
    if min_depth != 0:
        args.extend(['-d', str(min_depth)])
    
    if min_qual != 0.0:
        args.extend(['-q', str(min_qual)])
    
    # 布尔选项
    if exclude_missing:
        args.append('-e')
    
    if no_detailed:
        args.append('-D')
    
    if no_summary:
        args.append('-S')
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        vcf_stats_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\nVCF基因型统计分析被用户中断 | VCF genotype statistics analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"VCF基因型统计分析失败 | VCF genotype statistics analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv